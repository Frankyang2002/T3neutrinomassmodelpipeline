from pathlib import Path
import pytest
from Numerical.orchestration.ModelNames import model_label, shorten, long_name
from scripts.RenameT3ModelFolders import migrate


def test_classes():
    cases = {'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),'D':(3,1,2),'E':(3,3,2)}
    for label, dims in cases.items():
        for alpha in (-4,-2,-1,0,1,2):
            short = model_label(*dims,alpha)
            assert short == f'T3-{label}-' + ('m' if alpha<0 else 'p') + str(abs(alpha))
            assert shorten(long_name(short)) == short


def test_unknown_dimensions():
    with pytest.raises(ValueError):
        model_label(1,1,1,-1)


def test_migrate_dry_and_apply(tmp_path):
    folder = tmp_path/'Reports'/'output'/'models'/'T3_dS1_3_dS2_3_dF_2_alpha_p2'
    folder.mkdir(parents=True)
    (folder/'test.pdf').write_bytes(b'PDF data')
    assert migrate(tmp_path)['count'] == 1
    assert folder.is_dir()
    assert migrate(tmp_path,execute=True)['count'] == 1
    assert (folder.with_name('T3-E-p2')/'test.pdf').read_bytes() == b'PDF data'


def test_collision(tmp_path):
    root = tmp_path/'Reports'/'output'/'models'
    (root/'T3_dS1_3_dS2_3_dF_2_alpha_p2').mkdir(parents=True)
    (root/'T3-E-p2').mkdir()
    with pytest.raises(FileExistsError):
        migrate(tmp_path,execute=True)
