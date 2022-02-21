## LSS2 (Linear System Solvers)

LSS is C++ library written written in Visual Studio 2019 IDE (VS) on win-x64 platform using C++17 standard. It contains several linear system solvers plus some ODE and PDE solvers. It has been migrated from LSS using valarrays instead of vectors.
It also has some output functionalities.

### Requirements

LSS uses among others some CUDA solvers. Therefore NVidia's graphic device with CUDA library v11.3 (and above) is required.

### Dependency

* NVidia's CUDA library v11.3 (and newer)

### Setup

As the library is being developed and tested in VS I will describe the setup only for this IDE. For other IDEs the steps will be similar.

##### 0.Step
Download the library in zip or tar above and extract. 

##### 1.Step
Make sure to create CUDA project (>= v11.3) in VS

##### 2.Step
Open property pages for the newly created project and
1. under **Debugging** tab set **Environment** to point to CUDA binaries folder (in my case it is *PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.3\bin*)
2. under **VC++ Directories** tab set **Include Directories** to point to include folder of this library
3. under **VC++ Directories** tab set **Library Directories** to point to lib folder of this library
4. under **Linker->Input** tab add to **Additional Dependencies** *cusolver.lib*;*cusparse.lib*; *cublas.lib*; *cublasLt.lib*;*lss_debug.lib* (or *lss_release.lib* in case of Release configuration)
5. under **Linker->General** tab set **Additional Library Dependencies** to point to lib folder of this library
6. under **CUDA Linker->General** tab set **Additional Library Directories** to point to lib folder of this library

##### 3.Step
Place *lss_debug.dll*, *lss_debug.lib* (*lss_release.dll*, *lss_release.lib* in case of Release configuration) into your executable folder

##### 4.Step
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
If it does not compile or run ok make sure you did not miss any of the previous steps. Otherwise you can drop me a line on my github page.


### Contents

In this section I will try to describe how to work with the library on several examples. I will also list several conventions which a client has to follow in order to properly use the solvers inside the library. Some of these conventions will be relaxed as the developement progresses.

*  [Boundary and initial/terminal conditions conventions](./boundary-conventions.html)
*  [2nd order ODEs](./second-degree-ode.html)
*  [1D Heat type equations](./heat-type-1d.html)
*  [1D Wave type equations](./wave-type-1d.html)
*  [2D Heston type equations](./heston-type-2d.html)


### Support or Contact

Having trouble with LSS? You can open [issue](https://github.com/MichalSara99/lss/issues) or contact me on mail *michal.sara99@gmail.com* and I’ll try and help you sort it out.
