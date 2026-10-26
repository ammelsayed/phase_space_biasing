# phase_space_biasing

Study of the effects of **phase-space biasing** on V+jets production with
MadGraph5_aMC@NLO (MG5aMC): `pp → Z/W → leptons + 0,1,2 jets`, MLM merged,
13.6 TeV (`ebeam 6800`).

The idea: MG5aMC can bias the LO event generation so that more unweighted
events are produced in a region of phase space that is statistically
limited (e.g. the high-pT tail), while the physics content of the sample
is preserved by the corresponding event weights. This repository collects
the generation scripts and a set of biasing modules to quantify how much
statistics can be gained in the tails of the `ptl1`, `MET`, `HT` and
`ptj` distributions of V+jets samples.

## Repository layout

```
phase_space_biasing/
├── z012j_MLM.mg5                  # generation script, Z(→l+l-) + 0,1,2j (MLM)
├── w012j_MLM.mg5                  # generation script, W(→lν)   + 0,1,2j (MLM)
├── ptj_bias.f                     # verbatim copy of the shipped MG5aMC ptj_bias module
└── mg5_bias_modules/              # biasing modules (each one is a self-contained module dir)
    ├── ptj_bias/                  # the shipped module (reference/template)
    ├── ptl1_bias/                 # bias on the pt of the leading charged lepton
    ├── met_bias/                  # bias on the missing transverse energy
    ├── ht_bias/                   # bias on HT (scalar sum of visible pT)
    ├── ptl1_met_ht_bias/          # combined bias: ptl1 × MET × HT
    └── xsec_flat_bias/            # (experimental) bias on the cross section itself
```

Every module directory contains `<name>.f` (the Fortran implementation of
`bias_wgt`) and a `makefile`, following the layout of the shipped
`Template/LO/Source/BIAS/ptj_bias` module of MG5aMC.

## Generation scripts

`z012j_MLM.mg5` and `w012j_MLM.mg5` are self-contained MG5aMC scripts
(process definition + launch). Run them with:

```bash
mg5_aMC z012j_MLM.mg5
mg5_aMC w012j_MLM.mg5
```

Settings (as requested for this study):

| Setting | Value |
|---|---|
| Events | `set nevents 100_000` |
| Random seed | `set iseed 0` (chosen from the clock; no multi-seed scan) |
| Shower | `shower = OFF` (no Pythia8) |
| Detector | `detector = OFF` (no Delphes) |
| Analysis | `analysis = MadAnalysis5` |
| MadSpin / reweight | OFF |
| Merging | `ickkw 1` (MLM), `xqcut 30.0`, `maxjetflavor 5` |
| Cuts | `ptj 30`, `etaj 2.5`, `ptl1min 150`, `ptl2min 50` / `misset 50` (W), `etal 2.5` |
| PDF | `lhapdf`, `lhaid 331900` |

The process-definition block (`define`, `generate`, `output`, diagram
plotting) is uncommented, so each script regenerates the process from
scratch. Since Pythia8 and Delphes are OFF, the two
`./Cards_edited/pythia8_card_*.dat` / `delphes_card_*.tcl` lines are
commented out (they would otherwise be dangling answers to questions that
are no longer asked).

Prerequisites:

* A PDF set with id `331900` installed in LHAPDF.
* MadAnalysis5 for the analysis stage: install it (e.g. `install
  madanalysis5` from the MG5aMC prompt) and make sure `madanalysis5_path`
  is set in `mg5_configuration.txt`, otherwise `analysis = MadAnalysis5`
  cannot run.

## How MG5aMC phase-space biasing works

Documented in the [LOEventGenerationBias wiki
page](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias),
available since MG5aMC v2.5.0. Summary:

* The run_card entry `bias_module` selects the module (`None` by default,
  the name of a shipped module, or **an absolute path to a custom module
  directory**).
* Its double-precision parameters are set through the run_card entry
  `bias_parameters` (a python-style dictionary), or from the launch
  prompt with `set bias_parameter <name> <value>`.
* During unweighting, events are kept with a probability proportional to
  `|original_weight × bias_weight|`, so the bias reshapes where in phase
  space the unweighted events land.
* Each module declares, inside the Fortran code, the common block

  ```fortran
  double precision stored_bias_weight
  logical impact_xsec, requires_full_event_info
  common/bias/stored_bias_weight,impact_xsec,requires_full_event_info
  ```

  * `impact_xsec = .False.` : the bias is **only a sampling device**. The
    total cross section and the physical distributions are unchanged;
    the unweighted events come out *weighted* and each event's bias
    weight is recorded in the LHE file as `<rwgt><wgt id='bias'>`,
    so any analysis can undo the bias a posteriori. **All modules in
    this repository use this mode.**
  * `impact_xsec = .True.` : the bias is a genuine modification of the
    matrix-element weight (it changes the cross section and the
    distributions themselves).
  * `requires_full_event_info = .False.` : the module only needs the
    event kinematics `p(0:3,nexternal)` (faster). Set it to `.True.`
    (see the `PY8_CKKW` example on the wiki) if you need the full LHE
    event record, colors, helicities, etc.
* A custom module is a directory named after the module containing
  `<name>.f` (with the `subroutine bias_wgt(p, original_weight,
  bias_weight)` entry point) and a `makefile`. Parameters are declared
  with a comment inside the subroutine,

  ```fortran
  C  parameters = {'ptl1_bias_target_ptl1': 300.0,
  C               'ptl1_bias_enhancement_power': 4.0}
  ```

  MG5aMC parses it, creates the corresponding run_card entries and
  writes their values into `Source/BIAS/bias.inc`, which the module
  includes after the variable declarations.
* Cluster-mode caveat (from the wiki): a bias module must not read
  external files at runtime — compile everything in.

## The biasing modules

All modules have the functional form
`bias_wgt = (observable / target)^power` and leave events without the
relevant object unbiased (weight 1). Defaults below are what the run_card
gets if you don't override `bias_parameters`.

| Module | Observable | Default parameters |
|---|---|---|
| `ptj_bias` (shipped, copied here) | pt of the hardest jet | `target 1000`, `power 4` |
| `ptl1_bias` | pt of the hardest charged lepton (e±, μ±, τ±) | `ptl1_bias_target_ptl1 300`, `ptl1_bias_enhancement_power 4` |
| `met_bias` | MET = \|Σ pT of invisible particles (neutrinos)\| | `met_bias_target_met 300`, `met_bias_enhancement_power 4` |
| `ht_bias` | HT = scalar Σ pT of all visible particles | `ht_bias_target_ht 1000`, `ht_bias_enhancement_power 4` |
| `ptl1_met_ht_bias` | product of the three biases above | all six parameters above |
| `xsec_flat_bias` | (experimental) see below | `xsec_flat_bias_exponent 1.0`, `xsec_flat_bias_reference_weight 1.0` |

### xsec_flat_bias — biasing the cross section directly

Instead of biasing one kinematic distribution, this module tries to bias
the *cross section itself*: it uses the `original_weight` argument of
`bias_wgt`, which is proportional to the differential cross section at
that phase-space point, and inverts it:

```fortran
bias_wgt = xsec_flat_bias_reference_weight
           / |original_weight|^xsec_flat_bias_exponent
```

With exponent = 1 the combined weight `original_weight × bias_weight` is
flat, i.e. **every phase-space point is generated with the same rate** —
regions with a small cross section (which would normally be sampled
rarely) get as many events as the regions with a large one. `0 < exponent
< 1` gives a partial flattening, `0` recovers the unbiased generation.

It runs with `impact_xsec = .False.`, so the reported cross section and
the physical distributions are untouched; only the *population* of the
unweighted sample changes, and each event's bias weight is stored in the
LHE `<rwgt>` block for a-posteriori correction. Caveats (documented in
the module header): the flattening is in the generator's integration
measure (which already includes its own importance sampling), not
uniformly in a physical observable, and with exponent = 1 the bias
weights can be huge — if the unweighting efficiency degrades, lower the
exponent.

## Using the modules

Point `bias_module` at the absolute path of the module directory:

```bash
mg5_aMC> set bias_module /path/to/this/repo/mg5_bias_modules/ptl1_bias
mg5_aMC> set bias_parameters {'ptl1_bias_target_ptl1': 500.0, 'ptl1_bias_enhancement_power': 4.0}
```

or, per-parameter from the launch prompt:

```bash
mg5_aMC> set bias_parameter ptl1_bias_target_ptl1 500.0
```

To disable: `set bias_module None`.

## How to study phase-space biasing with this repository

1. **Reference sample**: generate `z012j_MLM` / `w012j_MLM` with
   `bias_module None`.
2. **Biased samples**: for each module, rerun the *same script* with the
   same `nevents` and cuts, changing only `bias_module` (and optionally
   its parameters).
3. **Compare**, per module:
   * the relative statistical uncertainty of the biased vs. reference
     distributions in the tails (`ptl1`, `MET`, `HT`, `ptj`), where the
     bias should buy you events;
   * the total cross section — with `impact_xsec = .False.` it must agree
     with the reference within uncertainties;
   * the central values — reweighting each biased event by
     `1/bias_weight` (or simply dividing the stored event weight by the
     `<wgt id='bias'>` entry) must reproduce the reference distributions,
     which is the basic closure test of the whole procedure.
4. **`xsec_flat_bias`**: additionally check how the unweighted events are
   spread over the phase space (population of normally rare regions) and
   the unweighting efficiency printed by MG5aMC, as a function of
   `xsec_flat_bias_exponent`.

## References

* [LOEventGenerationBias](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/LOEventGenerationBias)
  — official MG5aMC wiki documentation of the biasing feature, including
  the custom-module interface and the `PY8_CKKW` example used as the
  advanced reference.
* The shipped module `Template/LO/Source/BIAS/ptj_bias` of MG5aMC —
  copied verbatim in this repository (`ptj_bias.f` and
  `mg5_bias_modules/ptj_bias/`); every custom module here follows its
  structure.
* MG5aMC is available since v2.5.0 with this feature; the framework is
  described in J. Alwall *et al.*, "The automated computation of
  tree-level and next-to-leading order differential cross sections, and
  their matching to parton shower simulations", JHEP **07** (2014) 079,
  [arXiv:1405.0301](https://arxiv.org/abs/1405.0301).
* MadAnalysis5: E. Conte, B. Fuks, G. Serret, "MadAnalysis 5, a user
  framework for qualitative analysis of collider events", Comput. Phys.
  Commun. **184** (2013) 222, [arXiv:1206.1599](https://arxiv.org/abs/1206.1599).
