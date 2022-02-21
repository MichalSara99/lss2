# LSS2 (Linear System Solvers) - migrated from LSS
Linear System Solvers library. Written in C++17. It contains a few wrappers around CUDA cuSolver library functions plus some other well known solvers.
It also contains some ODE and PDE solvers.

## ODE Solver
* general 2nd degree two-point equation with variable coefficients (supports all Dirichlet, Neumann, Robin boundary conditions)
* All solvers support uniform and non-uniform grid with variable scaling


## PDE Solver
* 1D general heat equation with variable coefficients (in space and time dimension) (supports all Dirichlet, Neumann, Robin boundary conditions)
* 1D general wave equation with variable coefficients (in space and time dimension) (supports all Dirichlet, Neumann, Robin boundary conditions)
* 2D general Heston type model with variable coefficients (in space and time dimension)
 (support for Douglas-Rachford ADI, Craig-Sneyd ADI, Modified Craig-Sneyd ADI, Hundsdorfer-Verwer ADI)
* 3D general Heston-Hull-White type model with variable coefficients (in space and time dimension)
 (support for Douglas-Rachford ADI so far) - only experimental and not yet released in DLLs (to be released soon..)
* All solvers support uniform and non-uniform grid with variable scaling

## Requirement
* Library is being developed in VS2019 for win-x64 arch
* NVIDIA CUDA lib >= v11.3

## Usage
### Visual Studio:

#### 0.step
   Make sure to create CUDA project (>= v11.3) in VS
#### 1.step 
   Open property pages for the newly created project and
   1. under Debugging set Environment to point to CUDA binaries folder (in my case it is PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.3\bin\)
   2. under VC++ Directories set Include Directories to point to include folder of this library
   3. under VC++ Directories set Library Directories to point to lib folder of this library
   4. under Linker|Input add to Additional Dependencies cusolver.lib;cusparse.lib; cublas.lib; cublasLt.lib;lss2_debug.lib (or lss2_release.lib in case of Release configuration)
   5. under Linker|General set Additional Library Dependencies to point to lib folder of this library
   6. under CUDA Linker|General set Additional Library Directories to point to lib folder of this library
#### 2.step
   Place lss2_debug.dll,lss2_debug.lib (lss2_release.dll, lss2_release.lib in case of Release configuration) into your executable folder
#### 3.step
   Now you should be ready to use the library. Test it using following example


```cpp
// main.cu file 


#include <iostream>
#include <string>
#include <lss/utilities/range.hpp>
#include <lss/configs/discretization_config_1d.hpp>


int main()
{

    // number of space subdivisions:
    std::size_t const Sd = 100;
    // number of time subdivisions:
    std::size_t const Td = 100;
    // space range:
    auto const &space_range_ptr = lss::range_builder().lower(0.0).upper(20.0).build();
    // time range
    auto const &time_range_ptr = lss::range_builder().lower(0.0).upper(1.0).build();

    //// discretization config:
    auto const &discretization_ptr = lss::discretization_config_1d_builder()
                                        .space_range(space_range_ptr)
                                        .number_of_space_points(Sd)
                                        .time_range(time_range_ptr)
                                        .number_of_time_points(Td)
                                        .build();

    std::cout<<"Uniform space step value: "<< discretization_ptr->space_step()<<"\n";

    return 0;
}
```

## Output
* results can be put into stream or xml. Folder output also contains some python scripts that parse the resulting xmls and plot the solutions. 

## Some curves from ODE solver
Simple Two-point BVP (u''(x) = -2 with Dirichlet and Robin BC)
![Simple ODE equation](/outputs/pics/ode_bvp_neumann_robin.png)

## Some surfaces from PDE solver

Heat equation (Dirichlet BC) from implicit solver
![Pure heat equation](/outputs/pics/pure_heat_surf_cuda_euler_nonuniform_grid.png)

Wave equation (Dirichlet BC) from implicit solver
![Pure wave equation](/outputs/pics/pure_wave_surf_implicit_tlu_uniform_grid.png)

Wave equation (Neumann BC) from implicit solver
![Pure wave equation - neumann](/outputs/pics/pure_wave_surf_implicit_cuda_dev_nonuniform_grid.png)

Damped wave equation (Dirichlet BC) from implicit solver
![Damped wave equation](/outputs/pics/dumped_pure_wave_surf_implicit_dss_nonuniform_grid.png)

Advection equation from implicit solver
![Advection equation](/outputs/pics/advection_surf_implicit_tlu_nonuniform_grid.png)

Black-Scholes equation from implicit solver
![Black-Scholes equation](/outputs/pics/call_surf_bs_implicit_tlu_cn_uniform_grid.png)

Heston equation DR from implicit solver
![Heston equation DR](/outputs/pics/heston_tlu_dr_nonuniform_grid.png)

Heston Barrier equation DR from implicit solver
![Heston Barrier equation DR](/outputs/pics/heston_tlu_dr_uao_put_nonuniform_grid.png)

Heston equation MCS from implicit solver
![Heston equation MCS](/outputs/pics/heston_tlu_mcs_call_nonuniform_grid.png)

Heston equation HV from implicit solver
![Heston equation HV](/outputs/pics/heston_tlu_hw_call_nonuniform_grid.png)

SABR equation from implicit solver
![SABR equation](/outputs/pics/sabr_tlu_dr_call_nonuniform_grid.png)

Heston equation from explicit solver
![Heston equation expl](/outputs/pics/heston_explicit_euler_call_nonuniform_grid.png)

SABR dynamics
![SABR dynamics](/outputs/pics/sabr_tlu_dr_call_nonuniform_grid_stepping.png)
