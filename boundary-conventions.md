[Back](index.md)

## Boundary and initial/terminal conditions conventions

### Boundary Conditions

LSS supports three types of boundary conditions: **Dirichlet**, **Neumann** and **Robin**. All solvers except for specific ones (so far only Heston type) allow any combination of these. 

### ODE Solvers
In what follows ![solution_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B200%2C200%2C200%7D%20u%28x_%7B0%7D%29%20%7D) represents a solution of an ODE and ![boundary_point_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B200%2C200%2C200%7D%20x_%7B0%7D%20%7D) represents a point in which we evaluate a boundary condition. Let us consider following simple 2nd order ODE

![ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B200%2C200%2C200%7D%20%5Cfrac%7Bd%5E%7B2%7Du%28x%29%7D%7Bdx%5E%7B2%7D%7D%3D-2%20%7D)



#### Dirichlet Boundary
  
LSS assumes Dirichlet boundary of following form
  
 ![dirichlet_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28x_%7B0%7D%29%20%3D%20A%7D)

Specifically, if we have following Dirichlet boundary:

![dirichlat_exampl_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28x_%7B0%7D%29%20%3D%203.1415%7D)

then in code we write:

```cpp
  #include <lss/boundary/dirichlet_1d.hpp>

  auto const dirichlet =  3.1415; 
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::dirichlet_1d_builder().value(dirichlet).build();

```
 
  
#### Neumann Boundary

LSS assumes Neumann boundary of following form

 ![neumann_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7Bdu%28x_%7B0%7D%29%7D%7Bdx%7D%20&plus;%20A%20%3D%200%7D)
 
Specifically, if we have following Neumann boundary:

![neumann_exampl_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7Bdu%28x_%7B0%7D%29%7D%7Bdx%7D%20%3D%201%7D)

then in code we write:

```cpp
  #include <lss/boundary/neumann_1d.hpp>

  auto const neumann = -1.0;
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::neumann_1d_builder().value(neumann).build();

```

Note *-1.0* rather then *1.0*. We must pay attention to the aforementioned convention or we get completely different solution.

#### Robin Boundary

LSS assumes Robin boundary of following form

 ![robin_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7Bdu%28x_%7B0%7D%29%7D%7Bdx%7D%20&plus;%20Au%28x_%7B0%7D%29%20&plus;%20B%20%3D%200%7D)

Specifically, if we have following Robin boundary:

![robin_exampl_ode](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7Bdu%28x_%7B0%7D%29%7D%7Bdx%7D%20&plus;%202u%28x_%7B0%7D%29%20%3D%200%7D)

then in code we write:

```cpp
  #include <lss/boundary/robin_1d.hpp>

  auto const robin_first =  2.0; 
  auto const robin_second =  0.0;
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::robin_1d_builder()
                                  .values(robin_first,robin_second)
                                  .build();

```

### 1D PDE Solvers

In what follows ![solution_1d_pde](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28t%2Cx%29%7D) represents a solution of 1D PDE and ![boundary_point_1d_pde](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20x_%7B0%7D%7D) represents a point in which we evaluate a boundary condition and ![initial_point_1d_pde](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20t_%7B0%7D%7D) represents a time at which we evaluate an initial/terminal condition. Let us consider following 1D PDE (heat equation)

![pde_1d](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t%2Cx%29%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B%5Cpartial%5E%7B2%7Du%28t%2Cx%29%7D%7B%5Cpartial%20x%5E%7B2%7D%7D%7D)



#### Dirichlet Boundary
  
LSS assumes Dirichlet boundary of following form
  
 ![dirichlet_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28t%2Cx_%7B0%7D%29%20%3D%20A%28t%29%7D)

Specifically, if we have following Dirichlet boundary:

![dirichlet_exampl_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28t%2Cx_%7B0%7D%29%20%3D%202*t%7D)

then in code we write:

```cpp
  #include <lss/boundary/dirichlet_1d.hpp>

  auto const &dirichlet = [](double t) { return 2*t; }; 
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::dirichlet_1d_builder().value(dirichlet).build();

```
 
  
#### Neumann Boundary

LSS assumes Neumann boundary of following form

 ![neumann_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t%2Cx_%7B0%7D%29%7D%7B%5Cpartial%20x%7D%20&plus;%20A%28t%29%20%3D%200%7D)
 
Specifically, if we have following Neumann boundary:

![neumann_exampl_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t%2Cx_%7B0%7D%29%7D%7B%5Cpartial%20x%7D%20%3D%202t%7D)

then in code we write:

```cpp
  #include <lss/boundary/neumann_1d.hpp>

  auto const neumann = [](double t) { return -2.0*t; }; 
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::neumann_1d_builder().value(neumann).build();

```

Note *-2.0t* rather then *2.0t*. We must pay attention to the aforementioned convention or we get completely different solution.

#### Robin Boundary

LSS assumes Robin boundary of following form

 ![robin_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t%2Cx_%7B0%7D%29%7D%7B%5Cpartial%20x%7D%20&plus;%20A%28t%29u%28t%2Cx_%7B0%7D%29%20&plus;%20B%28t%29%20%3D%200%7D)

Specifically, if we have following Robin boundary:

![robin_exampl_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t%2Cx_%7B0%7D%29%7D%7B%5Cpartial%20x%7D%20&plus;%202tu%28t%2Cx_%7B0%7D%29%20%3D%200%7D)

then in code we write:

```cpp
  #include <lss/boundary/robin_1d.hpp>

  auto const robin_first =  [](double t) { return 2*t; }; 
  auto const robin_second = [](double t) { return 0.0; };
  
  // build boundary conditions from above function:
  auto const &boundary_ptr = lss::robin_1d_builder()
                                  .values(robin_first,robin_second)
                                  .build();

```

### Initial/Terminal Condition for 1D PDEs 

Initial/terminal condition follows similar pattern. 

![initial_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28t_%7B0%7D%2Cx%29%20%3D%20I%28x%29%7D)

which for

![initial_exampl_1d_pde](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20I%28x%29%20%3D%20max%280.0%2Cx-strike%29%7D)

translates into C++ code

```cpp
  #include <lss/pde/configs/heat_initial_data_config_1d.hpp>

  // strike
  auto const strike = 20.0;
  // terminal condition function (payoff of call):
  auto terminal_condition = [=](double x) { return std::max<double>(0.0, x - strike); };

  // build initial data config:
  auto const &heat_init_data_ptr =
        lss::heat_initial_data_config_1d_builder().condition(terminal_condition).build();


```

Note that builder *lss::heat_initial_data_config_1d_builder()* is used for either initial or terminal condition.

In case we have two initial/terminal conditions (as in case of wave type equations)

![initial_first_1d_pde_w](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20u%28t_%7B0%7D%2Cx%29%20%3D%20I%28x%29%7D)

![initial_second_1d_pde_w](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20%5Cfrac%7B%5Cpartial%20u%28t_%7B0%7D%2Cx%29%20%7D%7B%5Cpartial%20t%7D%20%3D%20J%28x%29%7D)

we will have for 

![initial_first_exampl_1d_pde_w](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20I%28x%29%20%3D%20%5Csin%7B%28%5Cpi%20*x%29%7D%7D)

![initial_second_exampl_1d_pde_w](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%7B%5Ccolor%5BRGB%5D%7B240%2C240%2C240%7D%20J%28x%29%20%3D%200.0%7D)

C++ code

```cpp
  #include <lss/pde/configs/wave_initial_data_config_1d.hpp>

  // initial condition:
  auto first_condition = [](double x) { return std::sin(pi() * x); };
  auto second_condition = [](double x) { return 0.0; };
  auto const wave_init_data_ptr = std::make_shared<wave_initial_data_config_1d>(initial_condition, second_condition);


```

[Back](index.md)


