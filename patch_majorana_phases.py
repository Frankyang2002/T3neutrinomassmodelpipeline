from pathlib import Path

path = Path("T3RGECorrectedYukawaPoint.py")
text = path.read_text(encoding="utf-8")

def replace_once(old: str, new: str) -> None:
    global text
    if old not in text:
        raise SystemExit("Patch target not found; no changes made. Missing:\n" + old[:200])
    text = text.replace(old, new, 1)

replace_once(
    "def required_high_scale_c5(\n    config: dict,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n) -> tuple[np.ndarray, np.ndarray]:\n",
    "def required_high_scale_c5(\n    config: dict,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n    alpha21: float = 0.0,\n    alpha31: float = 0.0,\n) -> tuple[np.ndarray, np.ndarray]:\n",
)

replace_once(
    "    target = build_neutrino_target(\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n        vev_gev=vev_gev,\n    )\n",
    "    target = build_neutrino_target(\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n        vev_gev=vev_gev,\n        alpha21=alpha21,\n        alpha31=alpha31,\n    )\n",
)

replace_once(
    "def build_rge_corrected_config(\n    c5_path: Path,\n    config: dict,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n) -> tuple[dict, dict]:\n",
    "def build_rge_corrected_config(\n    c5_path: Path,\n    config: dict,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n    alpha21: float = 0.0,\n    alpha31: float = 0.0,\n) -> tuple[dict, dict]:\n",
)

replace_once(
    "    low_target, high_required = required_high_scale_c5(\n        config,\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n    )\n",
    "    low_target, high_required = required_high_scale_c5(\n        config,\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n        alpha21=alpha21,\n        alpha31=alpha31,\n    )\n",
)

replace_once(
    "    updated = copy.deepcopy(config)\n    updated[\"ordering\"] = ordering.upper()\n    model = updated[\"t3\"]\n",
    "    updated = copy.deepcopy(config)\n    updated[\"ordering\"] = ordering.upper()\n    updated[\"alpha21\"] = float(alpha21)\n    updated[\"alpha31\"] = float(alpha31)\n    model = updated[\"t3\"]\n",
)

replace_once(
    "        \"Ordering\": ordering.upper(),\n        \"LightestMassEV\": lightest_mass_ev,\n",
    "        \"Ordering\": ordering.upper(),\n        \"LightestMassEV\": lightest_mass_ev,\n        \"Alpha21Rad\": float(alpha21),\n        \"Alpha31Rad\": float(alpha31),\n",
)

replace_once(
    "def write_rge_corrected_config(\n    c5_path: Path,\n    input_config_path: Path,\n    output_config_path: Path,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n) -> dict:\n",
    "def write_rge_corrected_config(\n    c5_path: Path,\n    input_config_path: Path,\n    output_config_path: Path,\n    *,\n    ordering: str = \"NO\",\n    lightest_mass_ev: float = 0.01,\n    alpha21: float = 0.0,\n    alpha31: float = 0.0,\n) -> dict:\n",
)

replace_once(
    "    updated, diagnostics = build_rge_corrected_config(\n        c5_path,\n        config,\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n    )\n",
    "    updated, diagnostics = build_rge_corrected_config(\n        c5_path,\n        config,\n        ordering=ordering,\n        lightest_mass_ev=lightest_mass_ev,\n        alpha21=alpha21,\n        alpha31=alpha31,\n    )\n",
)

replace_once(
    "    parser.add_argument(\n        \"--m-lightest\",\n        type=float,\n        default=0.01,\n    )\n",
    "    parser.add_argument(\n        \"--m-lightest\",\n        type=float,\n        default=0.01,\n    )\n    parser.add_argument(\n        \"--alpha21\",\n        type=float,\n        default=0.0,\n        help=\"Majorana phase alpha21 in radians.\",\n    )\n    parser.add_argument(\n        \"--alpha31\",\n        type=float,\n        default=0.0,\n        help=\"Majorana phase alpha31 in radians.\",\n    )\n",
)

replace_once(
    "    summary = write_rge_corrected_config(\n        args.c5,\n        args.config,\n        args.output,\n        ordering=args.ordering,\n        lightest_mass_ev=args.m_lightest,\n    )\n",
    "    summary = write_rge_corrected_config(\n        args.c5,\n        args.config,\n        args.output,\n        ordering=args.ordering,\n        lightest_mass_ev=args.m_lightest,\n        alpha21=args.alpha21,\n        alpha31=args.alpha31,\n    )\n",
)

replace_once(
    "    print(\"ordering =\", summary[\"Ordering\"])\n",
    "    print(\"ordering =\", summary[\"Ordering\"])\n    print(\"alpha21 [rad] =\", summary[\"Alpha21Rad\"])\n    print(\"alpha31 [rad] =\", summary[\"Alpha31Rad\"])\n",
)

path.write_text(text, encoding="utf-8")
print("Patched T3RGECorrectedYukawaPoint.py")
print("- alpha21 and alpha31 now propagate into the target")
print("- phases are stored in the fitted config and fit summary")
