"""Evaluate one-loop RGBeta InputForm on the canonical T3 numerical UV state."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import re
from typing import Any, Mapping
import numpy as np
from Numerical.core.BetaVector import beta_values_to_derivative,validate_beta_keys
from Numerical.core.State import SharedT3UVState,T3UVState
from Numerical.core.StateVector import StateVectorLayout,pack_uv_state

Number=float|complex
Value=Number|np.ndarray
_TOKEN_RE=re.compile(r'(?P<WS>\s+)|(?P<NUMBER>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)|(?P<NAME>[A-Za-z_$][A-Za-z0-9_$]*)|(?P<OP>[\+\-\*/\^\.,\(\)\[\]])')
@dataclass(frozen=True)
class _Token:
    kind: str
    text: str

def _tokenize(text: str) -> list[_Token]:
    tokens=[];position=0
    while position<len(text):
        m=_TOKEN_RE.match(text,position)
        if m is None:raise ValueError(f'Unsupported RGBeta InputForm syntax near {text[position:position+40]!r}.')
        position=m.end()
        if m.lastgroup!='WS': tokens.append(_Token(m.lastgroup or '',m.group()))
    tokens.append(_Token('EOF',''))
    return tokens

def _is_array(value: Value) -> bool:return isinstance(value,np.ndarray)

def _ordinary_multiply(left: Value,right: Value) -> Value:
    if _is_array(left) and _is_array(right):
        raise ValueError("RGBeta expression used '*' between two matrices. The current evaluator requires Matrix[...] or '.' for matrix multiplication.")
    return left*right

def _matrix_multiply(left: Value,right: Value) -> Value:
    if not (_is_array(left) and _is_array(right)):
        raise ValueError("RGBeta '.' product currently requires matrix operands.")
    return left@right

class _Parser:
    def __init__(self,text: str,environment: Mapping[str,Any]) -> None:
        self.tokens=_tokenize(text);self.index=0;self.environment=environment
    @property
    def current(self) -> _Token:return self.tokens[self.index]
    def advance(self) -> _Token:
        token=self.current;self.index+=1;return token
    def accept(self,text: str) -> bool:
        if self.current.text==text:self.advance();return True
        return False
    def expect(self,text: str) -> None:
        if not self.accept(text):raise ValueError(f'Expected {text!r}, found {self.current.text!r}.')
    def parse(self) -> Value:
        value=self.parse_expression()
        if self.current.kind!='EOF':raise ValueError(f'Unexpected trailing token {self.current.text!r}.')
        return value
    def parse_expression(self) -> Value:
        value=self.parse_term()
        while self.current.text in ('+','-'):
            op=self.advance().text;rhs=self.parse_term();value=value+rhs if op=='+' else value-rhs
        return value
    def parse_term(self) -> Value:
        value=self.parse_power()
        while self.current.text in ('*','/','.'):
            op=self.advance().text;rhs=self.parse_power()
            if op=='*':value=_ordinary_multiply(value,rhs)
            elif op=='/':
                if _is_array(rhs):raise ValueError('Division by a matrix is not supported.')
                value=value/rhs
            else:value=_matrix_multiply(value,rhs)
        return value
    def parse_power(self) -> Value:
        value=self.parse_unary()
        if self.accept('^'):
            exponent=self.parse_power()
            if _is_array(exponent):raise ValueError('Matrix exponent is not supported.')
            if _is_array(value):raise ValueError('Elementwise/matrix powers are not used by the current RGBeta T3 export.')
            value=value**exponent
        return value
    def parse_unary(self) -> Value:
        if self.accept('+'):return self.parse_unary()
        if self.accept('-'):return -self.parse_unary()
        return self.parse_primary()
    def parse_primary(self) -> Value:
        token=self.current
        if token.kind=='NUMBER':self.advance();return float(token.text)
        if token.text=='(':
            self.advance();v=self.parse_expression();self.expect(')');return v
        if token.kind=='NAME':
            name=self.advance().text
            if self.accept('['):value=self.call_function(name,self.parse_arguments_after_open_bracket())
            else:
                if name not in self.environment:raise ValueError(f'Unknown RGBeta symbol {name!r}.')
                value=self.environment[name]
            while self.accept('['):self.skip_index_arguments_after_open_bracket()
            return value
        raise ValueError(f'Unexpected RGBeta token {token.text!r}.')
    def parse_arguments_after_open_bracket(self) -> list[Value]:
        args=[]
        if self.accept(']'):return args
        while True:
            args.append(self.parse_expression())
            if self.accept(']'):break
            self.expect(',')
        return args
    def skip_index_arguments_after_open_bracket(self) -> None:
        depth=1
        while depth:
            token=self.advance()
            if token.kind=='EOF':raise ValueError('Unterminated RGBeta index bracket.')
            if token.text=='[':depth+=1
            elif token.text==']':depth-=1
    def call_function(self,name: str,args: list[Value]) -> Value:
        if name in ('Bar','Trans','Tr'):
            self._require_arity(name,args,1);v=args[0]
            if name=='Bar':return np.conjugate(v)
            if name=='Trans':return np.transpose(v) if _is_array(v) else v
            if not _is_array(v):raise ValueError('Tr[...] requires a matrix.')
            return np.trace(v)
        if name=='Matrix':
            if not args:raise ValueError('Matrix[...] requires at least one factor.')
            result=args[0]
            if not _is_array(result):raise ValueError('Matrix[...] factors must be matrices.')
            for factor in args[1:]:
                if not _is_array(factor):raise ValueError('Matrix[...] factors must be matrices.')
                result=result@factor
            return result
        raise ValueError(f'Unsupported RGBeta function head {name!r}.')
    @staticmethod
    def _require_arity(name: str,args: list[Value],count: int) -> None:
        if len(args)!=count:raise ValueError(f'{name} expects {count} argument(s), got {len(args)}.')

def evaluate_inputform_expression(expression: str,environment: Mapping[str,Any]) -> Value:
    return _Parser(str(expression),environment).parse()

def state_environment(state: T3UVState|SharedT3UVState) -> dict[str,Any]:
    state=state.validated()
    env={'gY':state.sm.gY,'g2':state.sm.g2,'g3':state.sm.g3,
         'yu':state.sm.yu,'yd':state.sm.yd,'ye':state.sm.ye,'lambdaH':state.sm.lambdaH}
    if isinstance(state,SharedT3UVState):
        env.update({'h':state.h,'MF':state.MF,'mSSq':state.mSSq,'lambdaS':state.lambdaS,
                    'lambda3':state.lambda3,'lambda4':state.lambda4,'lambda5':state.lambda5})
        return env
    for name in ('y1','y2','MF','mS1Sq','mS2Sq','lambdaS1','lambdaS2','lambdaH1',
                 'lambdaH2','lambda12','lambdaT3','lambdaH1Adj','lambdaH2Adj',
                 'lambdaS1Adj','lambdaS2Adj','lambda12Adj','lambda12Cross',
                 'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaS1bar2S2bar2',
                 'lambdaS1barS2S2bar2','lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross'):
        value=getattr(state,name)
        if value is not None:env[name]=value
    env.update(state.general_couplings)
    return env

def _validate_payload_metadata(payload: Mapping[str,Any],state: T3UVState|SharedT3UVState) -> None:
    metadata=payload.get('metadata',{})
    if not isinstance(metadata,Mapping):raise ValueError('RGBeta payload metadata must be an object.')
    rep=state.representation.validated()
    expected={'SharedScalar':rep.shared_scalar,'dS1':rep.d_s1,'dS2':rep.d_s2,'dF':rep.d_f,'alpha':rep.alpha}
    actual={'SharedScalar':bool(metadata.get('SharedScalar',False)),'dS1':metadata.get('dS1'),
            'dS2':metadata.get('dS2'),'dF':metadata.get('dF'),'alpha':metadata.get('alpha')}
    differences=[f'{k}: payload={actual[k]!r}, state={v!r}' for k,v in expected.items() if actual[k]!=v]
    if differences:raise ValueError('RGBeta payload does not match the numerical state: '+'; '.join(differences))
    convention=metadata.get('ReportBetaConvention')
    required='16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'
    if convention!=required:raise ValueError(f'Unsupported RGBeta report-beta convention. Expected {required!r}, got {convention!r}.')

def evaluate_rgbeta_payload(payload: Mapping[str,Any],state: T3UVState|SharedT3UVState) -> dict[str,Value]:
    status = payload.get('status')
    fmt = payload.get('format')
    general = (fmt == 'general-t3-quadratic-rgbeta-experimental-v1')
    if general:
        if status != 'SuccessCandidateMassBeta':
            raise ValueError('General RGBeta payload must have SuccessCandidateMassBeta status.')
        metadata = payload.get('metadata')
        if not isinstance(metadata, Mapping) or metadata.get('GeneralScalarStage') != 'UV':
            raise ValueError('UV evaluator requires a general UV-stage RGBeta payload.')
        if not isinstance(state, T3UVState):
            raise ValueError('General ordinary-T3 payload requires T3UVState.')
    elif status != 'Success':
        raise ValueError('RGBeta payload status is not Success.')
    state=state.validated();_validate_payload_metadata(payload,state)
    betas=payload.get('report_betas')
    if not isinstance(betas,Mapping):raise ValueError('RGBeta payload does not contain a report_betas object.')
    _, layout = pack_uv_state(state)
    from Numerical.core.BetaVector import validate_beta_keys
    validate_beta_keys({name: None for name in betas}, layout)
    env=state_environment(state)
    evaluated={str(name):evaluate_inputform_expression(str(expr),env) for name,expr in betas.items()}
    _,layout=pack_uv_state(state)
    validate_beta_keys(evaluated,layout)
    return evaluated

def load_rgbeta_payload(path: Path) -> dict[str,Any]:
    payload=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(payload,dict):raise ValueError('RGBeta JSON top level must be an object.')
    return payload

def derivative_from_rgbeta_payload(payload: Mapping[str,Any],state: T3UVState|SharedT3UVState) -> tuple[np.ndarray,StateVectorLayout]:
    state=state.validated();_,layout=pack_uv_state(state)
    evaluated=evaluate_rgbeta_payload(payload,state)
    return beta_values_to_derivative(evaluated,layout,divide_by_loop_factor=True),layout
