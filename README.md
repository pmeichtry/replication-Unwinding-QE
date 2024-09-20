# Unwinding Quantitative Easing: Replication package

Replication files for the main results of Cantore, C. & Meichtry, P. (2024). Unwinding quantitative easing: State dependency and household heterogeneity. _European Economic Review, 170_, 104865.\
[[Published article](https://www.sciencedirect.com/science/article/abs/pii/S0014292124001946)] &nbsp; [[Complete version (ungated)](https://pmeichtry.github.io/Papers/UnwindingQE_CantoreMeichtry_paper.pdf)] &nbsp; [[Online Appendix](https://pmeichtry.github.io/Papers/UnwindingQE_CantoreMeichtry_OnlineAppendix.pdf)]

MATLAB and Dynare codes for New Keynesian model with borrowers and savers, short-term and long-term bonds, and a zero lower bound on nominal interest rates.

For any questions or comments, please contact pascal.meichtry@banque-france.fr.

# Structure
**codes**: codes to run the model simulations and create figures and tables
* **main_replication.m**: master file to replicate the main results
* **\functions**: auxiliary functions
* **\dynareOBC**: dynareOBC toolbox ([dynareOBC GitHub repository](https://github.com/tholden/dynareOBC/releases))

**simul_results**: Stored results from model simulation

**outputs**: Stored figures and table


See **[README.pdf](README.pdf)** contained in the package for more details on the replication.
