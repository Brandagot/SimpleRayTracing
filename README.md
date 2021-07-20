    ---
Original Author: Dr Franck P. Vidal <br />
Institute: School of Computer Science and Electronic Engineering, Bangor University, UK <br />
Title: ICE4131 -- High Performance Computing <br />
Subtitle: A Simple ray tracer to parallelise using PThread, OpenMP, MPI and CUDA. <br />
---
<br />

# SimpleRayTracing

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2e8c4fd913234f2d880c43716c17cea9)](https://app.codacy.com/manual/effepivi/SimpleRayTracing?utm_source=github.com&utm_medium=referral&utm_content=effepivi/SimpleRayTracing&utm_campaign=Badge_Grade_Dashboard)

Originally used as an example of application for Dr. Franck P. Vidal's HPC module. The code here has been forked to be modified as a software implementation of X-ray attenuation on the CPU as a part of my MSc thesis.

## Test Log

This section demonstrates some images created during testing.

At the moment, I have made some modifications to remove the dragon and changed the way distance is used to set the shade. The x-ray source is now in the same position as the source.

![](out/x-ray-10-threads-1512x1512-no_dragon.jpg)

Progress made, now using L-buffer to produce the following image.

![](out/attenuation_model_2048x2048_lbuffer.jpeg)
