# Purpose of this repository
This repository will provide material that has been omitted from a written document on a math project (such as .gif files and code).

# Short excursion into the history of MCMC

\section{Appendix}
In the appendix, I would like to make a small excursion into the history of MCMC methods.
The purpose of this is to provide context for the implementation of MCMC methods and to try to appreciate how far its come since its first applications in statistical mechanics.

\subsection{Supplementary material}
Having origins in statistical mechanics, we see how MCMC methods have evolved over time to become a general tool to solve many optimization problems of interest. We will first explore the work of Metropolis et. al in their application of MCMC methods. Their choice of procedure can be a little ambiguous for a reader unacquainted with notions of equilibrium in physical systems, so for clarity's sake, the work of Kirkpatrick et. al will be explored to understand the relation between statistical mechanics and optimization problems. 


\section{History of MCMC Methods}
\subsection{Metropolis et. al (1953) \hyperref[sec:subref]{$[2]$}}
\label{sec:metropolis}
Metropolis et. al published a paper titled "Equation of State Calculations by Fast Computing Machines" in 1953 in the \textit{Journal of Chemical Physics}. The quantity of interest to be evaluated was defined as 

\begin{equation*}
    \bar{F} = \frac{\int F(\theta)exp\{-E(\theta)/kT\}d\theta}{\int       exp\{-E(\theta)/kT\}d\theta}
\end{equation*}

\noindent Heuristically, the quantity of interest is an equilibrium state. The integral is computed on $\mathbb{R}^{2N}$. It is evaluated over $\theta$, the set of N particles on $\mathbb{R}^2$. E is denoted as the energy and is defined by the potential function $V$ and the Euclidean distance $d_{ij}$ between particles $i, j \in \theta$.

\begin{equation*}
    E(\theta) = \frac{1}{2} \sum_{i=1}^{N} \sum_{\substack{j=1 \\ j\neq i}}^{N} V(d_{ij})\\
\end{equation*}

\noindent In the quantity of interest is the Boltzmann distribution $exp\{-E(\theta)/kT\}$. Later, we will see that many of the invariant distributions in application have a similar flavor to the Boltzmann distribution.\setlength{\parskip}{6pt}

\noindent For now, define $Z(T) = \int exp\{-E(\theta)/kT\}d\theta$ to be a normalizing constant of $\bar{F}$. 
The integral $Z(T)$ is usually an analytically intractable problem. Conventional Monte Carlo methods  which relied on random sampling are inadequate for the purpose of evaluating $Z(T)$, as it is difficult to sample from large dimensional problems.\setlength{\parskip}{6pt}

\noindent Metropolis et. al proposed a random walk mechanism which caused slight perturbations on the particles of the system. The perturbations for particle $i$ with position $(x_i, y_i)$ on a usual 2D-grid are given as
\begin{equation*}
    \begin{aligned}
        x_{i}^{proposal} & = x_{i} + \sigma\xi_{1i}\\
        y_{i}^{proposal} & = y_{i} + \sigma\xi_{2i}
    \end{aligned}
\end{equation*}

\noindent where $\sigma$ is the maximum allowed perturbation with adjustments in direction and magnitude by $\xi_{1i}$, $\xi_{2i}$ $\sim$  $Unif(-1, 1)$. The proposed movement would put the particle within a square of side $2\sigma$ from its original location. \setlength{\parskip}{6pt}

\noindent The energy difference is computed between the proposed configuration $c^*$ and the current configuration $c$; the proposal is accepted with probability

\begin{equation*}
    \alpha(c, c^*) = \text{min}\{1, exp(-\Delta E/kT\}
\end{equation*}

\noindent According to probability $\alpha(c, c^*)$, if $\Delta E \leq 0$, then the proposal is accepted as the equilibrium configuration of particles tend toward lower energy states. If $\Delta E \geq 0$, then the movement is accepted with probability $exp\{-\Delta E/kT\}$. Heuristically, the procedure mimics the behavior of a system of particles in a heat bath with temperature T.

\noindent If the proposal is rejected, then the current configuration is replicated in the random walk. The procedure is done by proposing a movement for random particles in the system for $M$ iterations. Under conditions of the Markov Chain being irreducible and reversible, the algorithm produces a result that approximates the desired integral by taking the average of the property of interest after each proposal

\begin{equation*}
    \hat{\bar{F}} \hat{=} \frac{1}{M} \sum_{j=1}^M F_j
\end{equation*}

\subsection{Kirkpatrick, Gelatt, Vecchi (1983)
\hyperref[sec:subref]{$[3, 5]$}}
\label{sec:kirk}

In a paper titled "Optimization by Simulated Annealing", the authors describe the relationship between statistical mechanics and combinatorial optimization and how the Metropolis Algorithm readily translates over to solving optimization problems.

In many statistical mechanic problems, there is in order of size of $10^{23}$ atoms per cubic centimeter in fluids and solids. Only small samples of these atoms can be observed in an experimental setting and their behavior must be taken as probabilistic. Particle positions in a configuration are weighted by the Boltzmann distribution, as seen similarly in Metropolis et. al (1953). Tinkering with the temperature T, i.e. simulated annealing, one could attempt to observe phase transitions of particles from gas, fluids, and solids. An interest of statistical mechanics is to uncover the behavior of atoms in the limit of lowest temperatures, which confers lowest energy configurations.

The process of simulated annealing readily translates to optimization. Instead of thinking about the behavior of particles in their lowest energy configurations, we think of the behavior of a Markov Chain in the setting of minimizing a cost function. The process of optimization is akin to heating the system to an unstable state with high temperature, then slowly decreasing the temperature and observing the steady state of each configuration of the system until the steady state freezes, i.e. becomes invariant. In accord to the analog to simulated annealing, the system that is being optimized eventually reaches an invariant, near-optimal position with minimized costs. The number of re-arrangements of the system at each temperature is called an annealing schedule.

An example of simulated annealing applied to optimization is the traveling salesman problem (TSP).
The problem lies in optimizing the cost of traveling between N cities. We believe that cost of traveling is proportional to distance travelled, so minimizing costs is proportional to minimizing distance traveled. 

Essentially, the positions of each city on a map are recorded into a list. The salesman chooses an initial tour from the list, wherein he starts at some city at the beginning of the tour, visits every other city once, and returns to the initial city. The initial tour is likely far from optimal; the tour is perturbed in some manner and the perturbation is accepted or rejected according to a probability that follows the form of $exp(-\beta\ell(x)$) (note that this has a similar flavor to the Boltzmann distribution), where $\beta$ is an interaction parameter inversely related to temperature and $\ell(x)$ is the length of a tour $x$. By starting with a low $\beta$ (high temperature), we slowly increase $\beta$ (decrease temperature), and then we follow the procedure of simulated annealing until the tour can no longer be perturbed in a way that decreases the distance traversed any farther, i.e. the cost function is eventually minimized. The choice of perturbation may vary; one example may be to reverse some portion of the tour between two cities on a given tour or some other type of permutation.



\subsection{Putting the Hastings in Metropolis (1970)
\hyperref[sec:subref]{$[4]$}}
\label{sec:hastings}

In light of a more robust procedure that handles high dimensional problems with respect to optimization, many problems in statistics became computationally approachable. In "Monte Carlo Sampling Methods Using Markov Chains and Their Applications", Hastings used the MH-algorithm to work through problems that were previously considered intractable. Hastings' work demonstrated the broad applicability of MCMC methods beyond statistical physics.

\section{References}
\label{sec:References}

\noindent \href{www.nyu.edu/classes/tuckerman/stat.mech/lectures/lecture_26/node2.html}{[1]} “Exact Solutions of the Ising Model in 1 and 2 Dimensions.” NYU. \hyperref[sec:exact]{$\hookrightarrow$}

\noindent \href{https://pages.uoregon.edu/dlevin/MARKOV/markovmixing.pdf}{[2]} "Markov Chain Monte Carlo: Metropolis and Glauber Chains." Markov Chains and Mixing Times, by David Asher Levin et al., American Mathematical Society, 2017, pp. 37-45. 

\noindent \href{http://itf.fys.kuleuven.be/~enrico/Teaching/monte_carlo_2014.pdf}{[3]} Carlon, Enrico. “Computational Physics: Advanced Monte Carlo Methods.” Institute for Theoretical Physics. \hyperref[sec:Wolff Algorithm]{$\hookrightarrow$}

\noindent \href{https://services.math.duke.edu/~rtd/EOSP/EOSP2E.pdf}{[4]} Durrett, Richard. "The Metropolis-Hastings Algorithm." Essentials of Stochastic Processes. pp. 35-38.

\noindent [5] Based on lecture notes from MATH 5652  taught by Professor Arnab Sen in Fall 2017 semester. \hyperref[sec:sen]{$\hookrightarrow$}

\subsection{References for the Appendix}
\label{sec:subref}

\href{https://arxiv.org/pdf/0808.2902.pdf}{[1]} Robert, Christian, and George Casella. “A Short History of Markov Chain Monte Carlo: Subjective Recollections from Incomplete Data.” Statistical Science, vol. 26, no. 1, 2011, pp. 102–115. %JSTOR, JSTOR, www.jstor.org/stable/23059158.
% Outlines history of MCMC methods; exploring this below

\noindent \href{http://bayes.wustl.edu/Manual/EquationOfState.pdf}{[2]} Metropolis, N., Rosenbluth, A., Rosenbluth, M.,
Teller, A. and Teller, E. (1953). Equations of state calculations by fast computing machines. J. Chem. Phys. 21 1087–1092. \hyperref[sec:metropolis]{$\hookrightarrow$}
% Earl(iest?) demonstration of MCMC method in statistical mechanics

\noindent \href{http://www.jstor.org.ezp2.lib.umn.edu/stable/1690046?seq=1#page_scan_tab_contents}{[3]} Kirkpatrick, S., et al. “Optimization by Simulated Annealing.” Science, vol. 220, no. 4598, 1983, pp. 671–680. %JSTOR, JSTOR, www.jstor.org/stable/1690046.
% Used to show statistical mechanics -> Optimization
\hyperref[sec:kirk]{$\hookrightarrow$}

\noindent \href{https://pdfs.semanticscholar.org/a430/cb602f7d47c21cc4fa94a351ec0c4a9f1fbd.pdf}{[4]} Hastings, W. K. “Monte Carlo Sampling Methods Using Markov Chains and Their Applications.” Biometrika, vol. 57, no. 1, 1970, pp. 97–109. %JSTOR, JSTOR, www.jstor.org/stable/2334940.
% Showed that MCMC could be used to target statistical distributions
\hyperref[sec:hastings]{$\hookrightarrow$}

\noindent \href{https://services.math.duke.edu/~rtd/EOSP/EOSP2E.pdf}{[5]} Durrett, Richard. "The Metropolis-Hastings Algorithm." Essentials of Stochastic Processes. pp. 35-38.
% Used to identify Boltzmann Distribution probability weight for TSP
\hyperref[sec:kirk]{$\hookrightarrow$} 


# Monte Carlo Markov Chain - Ising Model
The Ising Model is a model that tries to capture the properties of interacting particle systems, namely in ferromagnetic materials. The model allows us to study topics such as energy, magnetization, and phase transitions. The Metropolis Algorithm (Glauber Dynamics) is an algorithm that simulates the Ising Model in thermodynamic equilibrium, and allows us to study properties of interest within the context of a computer simulation.

# Clustering Algorithms - Wolff Algorithm

The Wolff Algorithm was introduced by Uli Wolff in his 1989 paper titled "Collective Monte Carlo Updating for Spin Systems", which built upon the Swedsens-Wang non-local clustering algorithm. We would like to implement it as it overcomes some of the difficulties facing the single-site algorithm. 
