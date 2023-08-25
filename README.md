# GlaDS simulations and analysis for laminar--turbulent flux parameterization

## Structure

### 'analysis/'
Top-level analysis

### `glads/`:

`data/`: Model input data

`get_para.m`: Global default parameter settings

**Main cases**:

`00_synth_forcing/`

`01_kan_forcing/`

**Additional cases**:

`00a_shmip_forcing`: Melt forcing identical to SHMIP case D3

`01a_kan_adj_forcing`: KAN forcing with melt volume equal to SHMIP D3

`01b_kan_forcing_ev`: Englacial storage parameter 1e-5 (default: 1e-4)

`01c_kan_diurnal`: KAN forcing with prescribed diurnal melt-rate variations

`XXx_..._trough` and `XXx_..._valley`: Trough and valley bed topography

**Supplementary experiments**:

`S00_mesh_refinement`: Steady-state mesh refinement test

`S01_parameter_sensitivity`: Parameter sensitivity experiment

