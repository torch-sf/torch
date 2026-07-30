# Developer Guide

The Torch code is currently available on [GitHub](https://github.com/torch-sf/torch).


## Massive Star Cluster Simulations

Clouds with initial mass $\geq 10^6~\rm M_\odot$ are considered to be massive, and simulating the clusters that form from them requires special considerations. Here we provide a brief guide on how to simulate massive cluster formation in `Torch`.

### Mass Agglomeration

In massive clouds, the star formation efficiency will be between $80-100%$. This gives you an estimate of how many stars your cloud will form. A $10^6~\rm M_\odot$ cluster with a fully sampled IMF will easily produce over 1 million stars. It is highly recommended to use **mass agglomeration** in massive cluster simulations to reduce the computational load on the dynamics code. Torch has been run before with 800k stars, but it is very slow. Consider your science goals, how far you want to run your simulation, and what you want to study when deciding what your agglomeration mass should be. For example, if you want to get to gas expulsion, agglomerate stars such that your cluster will contain no more than ~250k particles. Another thing to consider is your feedback mass. If you choose an agglomeration mass of 4, you can have agglomerated stars with masses up to 8, so you must set your feedback mass above 8. 

>We recommend an agglomeration mass of $4\rm~M_\odot$, which will reduce the number of stars produced by $\sim80%$. 

### Feedback Mass

In massive cluster formation, there are so many massive stars that the feedback from low mass stars becomes negligible. Therefore, using a higher feedback mass cutoff reduces the computational load while only losing a small fraction of the total feedback energy budget. This is particularly important when using `FERVENT` for radiative feedbcak, as the time complexity scales as $O(N_\star^2)$.
>We recommend using a feedback mass of $20~\rm M_\odot$, and using `VETTAM` over `FERVENT`.

### Resolution

It can be tempting to start a high resolution run, but when massive stars flood the grid, your simulation will start to grind to a halt. 
>We recommend a smallest cell width no less than $\Delta x\sim0.3~\rm pc$ to ensure the run progresses steadily when $>10^6\rm~K$ gas forms from stellar winds. 

### Mass Loading

There will be a lot of feedback stars, and therefore a lot of stellar winds. The temperature of the gas sets the timestep of your simulation, so it is important to limit the temperature of stellar wind bubbles with mass loading. 
> We recommend a mass loading target temperature of $10^6\rm~K$. You can consider going lower, but be mindful that the peak of the cooling curve is at $\sim 3\times10^5\rm~K$, so your winds may cool artificially fast. 

### Stellar Evolution

`SeBa` in `Torch` can take up to 30 minutes to complete a single evolution step when ran with >200k stars. This is an ongoing issue that is being worked on. For now, the option to use interpolated values from a pre-tabulated `SeBa` steps is available. This is particularly useful for massive cluster formation, when simulations typically only run for 1 Myr or less. Stellar properties hardly change on this timescale.
> [OPTIONAL] If SeBa is your bottleneck, consider using static stellar evolution, available on the `feature/static-se` branch. 

### Initial Conditions

Very massive and dense clouds cannot be initialized in hydrostatic equilibrium. It is likely a pressure shockwave will form if your cloud isn't embedded in a high enough density medium. Extending the turbulent velocities to the entire box can also help with this. 
> We recommend using a high enough density ($100\rm~cm^{-3}$) ambient medium and applying the turbulent velocities to the entire domain if your initial cloud is very dense and unstable. 