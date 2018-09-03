**Abstract**

Milonga is a free computer code that solves the steady-state multigroup
neutron diffusion equation using either a finite-volumes or a
finite-differences scheme. Not only is it designed to cope with common
reactor geometries but also to parametrically study the effect of one or
more parameters in order to optimize some aspect of the reactor design.
The code is especially designed emphasizing flexibility in the way the
geometry and the cross-sections distributions are entered, including
dependence on arbitrary parameters such as temperatures, burn-up,
poisons, etc. This information can be entered as algebraic expressions,
multidimensional interpolated tables or from data given by external
codes trough shared memory objects. Milonga can handle a wide variety of
one, two and three-dimensional cases, from simple idealized problems up
to common reactor configurations, including xenon effects and coupled
calculations with thermal-hydraulics and control-logic codes. For
example, milonga can be used to solve a case with analytical solution
using its built-in algebraic and differential functions, and then easily
compare it to the discrete solution along with the corresponding CPU
times as a function of the number of spatial cells. The code output is
completely defined by the user through the input. The eigenvalue problem
is solved by the SLEPc+PETSc libraries, so the computational
implementation can be scaled virtually up to the limit of current
available hardware, and hopefully for many years to come. Moreover, not
only can milonga utilize user-provided *ad-hoc* numerical routines, but
it is also released under the GNU Public License so further scalability
and improvements may be introduced at will.




![image](sombrero){width="4cm"}


Preface to version 0.1 {#preface-to-version-0.1 .unnumbered .unnumbered}
======================

This document describes version 0.1 of the free nuclear reactor analysis
code milonga that is part of the developments made under my PhD thesis
in Nuclear Engineering. In particular, this version of the code and
especially of this document the very first not-so-public release and is
aimed at getting as much feedback as possible. The objective is to
improve the code to get into a version 0.9 that hopefully should
converge into a fully-usable version 1.0 in a nearby future. By
not-so-public I mean that this document was distributed mainly between
friends and colleagues. Thus, if you have received this document
directly from me---meaning you are a friend of mine---then I firmly ask
you to say something back. On the other hand, if you have got this
document from other indirect source, you are kindly encouraged to give
your comments back. These may be a linear combination of suggestions,
corrections, bug reports and improvements, in increasing order of
desirability.


This version of the code is not entirely complete in the sense that
there are some kind of problems it cannot solve in its present state.
Particularly, it does not fully support three dimensional calculations,
mainly because of its lack of management of boundary conditions.
Chapter [4](#cap:examples){reference-type="ref"
reference="cap:examples"} gives a wide variety of examples that this
version can solve. Completion of these capabilities is expected for the
next version.

In the same sense, this version of the documentation is neither entirely
complete. There are some sections marked as "to be done" meaning the
code is lacking that particular feature. There are also some sections
containing the text "to be explained". In this case, the code is able to
perform the task but it is not documented. This is closely related to
the fact that milonga is part of a set of codes that share a common
framework, conforming a suite of engineering codes called wasora that is
also part of my PhD thesis and of course are under development.


I would like to point out here in this preface that milonga is based on
the PETSc library to handle the matrices that are the base of the
multigroup diffusion problem formulation. This library provides very
efficient methods for creating, accessing and operating matrices. It is
especially designed to work in parallel using MPI. The current version
0.1 of milonga does not take any advantage of the parallelization
facilities provided by PETSc. In fact, this first version is very
computationally demanding because efficiency was not the main focus.
Optimization is expected to drive my attention as of version 1.0.


Finally, I want to say that milonga---and indeed the whole wasora
suite---was conceived as free software from scratch. Not only are
hacking and distribution under the terms of the GNU General Public
License allowed, but also encouraged. And of course, the kind gesture of
giving feedback to the original author will be highly appreciated.


Buenos Aires, July 2011



Introduction {#cap:introduction}
============



An engineer has to think without doing anything fifty percent of the
time,\
and do without thinking the other fifty percent of the time.\
And a good engineer knows when to think without doing\
and when to do without thinking.


*Fabián. J. Bonetto, PhD in Nuclear Engineering\
comment to the author during an engineering thesis advisory, 2007*


The first thing that should be said in relation to this code is that it
is an academic project part of a PhD thesis in Nuclear Engineering. This
thesis is not about software development, but about its usage as a
design optimization tool. That is to say, milonga is not the goal but
the mean. The second statement is that as such, it is a personal
research project written from scratch by someone who works in the
nuclear industry and is not fully dedicated to the academic environment.

Having said these two things, a corollary in the form of a warning is
issued: you should not expect a lot from milonga. In particular, it
should not be compared to commercial reactor codes, nor directly applied
to real cases without a deep comprehension of both its design philosophy
and the mathematics behind it.

No propaganda about the convenience of using milonga instead of other
codes will be given in this document because of three reasons. Firstly,
if this code does provide some interesting features that maybe other
program does not, its drawbacks far overcome the benefits. Secondly,
there is no need to do marketing because there is no commercial interest
behind milonga, as it is free software both in the speech and in the
beer senses. The implications are explained in
sections [1.2](#sec:basis){reference-type="ref" reference="sec:basis"}
and [1.3](#sec:lincense){reference-type="ref" reference="sec:lincense"}.
And last but not least, there is a problem with the conjunction used in
the first sentence of this paragraph. As in many other aspects, it is
not milonga *or* other codes, but milonga *and* other codes.

The wasora suite
----------------

Milonga is part of wasora, which stands for Wanna-be Advanced Suite for
Optimization and Reactor Analysis, that in turn is part of a PhD thesis
work. It consists of a general code framework that provides routines and
methods for engineering programs that share a common basis and may work
coupled. In particular it provides input-file keyword parsing routines
and an efficient method for parsing and evaluating algebraic
expressions. It also includes routines for accessing standard
mathematical functions, interpolation, root finding, integration and
differentiation. Most of the mathematical methods provided by wasora are
implemented by the free GSL library [@gsl-manual]. Access to shared
memory objects and semaphores for data exchange and coupling is also
provided, along with text-based output through files and error handling
routines. From the user's point of view, this just means that codes
built on top of wasora work in a similar way---for example their inputs
look all alike---and that they share certain common functionality.

The whole suite is under development but, as the current version of
milonga shows, the framework is usable. As of July 2011, the engineering
codes that are either under development or under consideration are:

colach

:   a control logic analysis code to quickly design, implement and test
    control algorithms in the time domain. Besides regular linear
    control function such as lags and integrators, matrix-vector
    operations, arbitrary transfer functions and fuzzy logic rules can
    be entered. A real-time version capable to operate with data
    acquisition hardware using COMEDI is under consideration.

mochin

:   a dynamical systems solver code to numerically integrate sets of
    differential-algebraic initial-value equations. It uses library IDA
    from the SUNDIALS suite as a back-end for time-advancing the set
    DAEs and the wasora framework for parsing the equations and probably
    coupling the code to other processes. Real-time simulation is also
    possible.

besssugo

:   a graphical-interface code to generate graphical time-dependent
    representations of the results obtained with engineering codes,
    either part of the wasora suite or not. Output may be real-time
    graphics in a screen or a series of frames to build a video
    afterward.

prime

:   a non-linear numerical optimization code. Under consideration.

cingi

:   a multi-point reactor kinetics model solver. Under consideration.

Further information and updates be obtained by accessing

<https://www.seamplex.com/wasora>


Of course, milonga---whose description is this tiny document---is also
part of the suite. A great deal of milonga's functionality directly
depends on the routines and methods provided by wasora, so its
documentation should complement this description. Sadly, at the time of
writing the present document for the first version of milonga, there is
no consistent documentation for the wasora framework available. Some
explanation about the basic usage---from the user's point of view---is
given in this document, but probably there may be some gaps either in
the information given here or when trying to understand some
particularities of the code. Hopefully, future releases will contain a
complete set of documentation covering every aspect of both wasora and
milonga.

Milonga design basis {#sec:basis}
--------------------

The quote of the beginning of this chapter was told privately to me by
Dr. F. J. Bonetto while he was my BE thesis advisor [@theler2007]. Later
on, he advised my Master's thesis [@theler2008] and he is now[^1] one of
my PhD thesis advisors. While years passed by, I have been able to
understand the rationale behind this statement. It is important to note
that the implementation of this phrase actually gave rise to remarkable
results all in all. Indeed, the development of milonga is an actual
consequence of this advise. This project has a lot of time spent just in
thinking rather than programming and, of course, the other way round.

Many of milonga's features were coded just after wondering what features
I wanted commercial programs to have and how I would have liked to work
with canned codes available in the nuclear industry. Indeed, one of the
reasons of the high component of "thinking without doing" comes from the
fact that working all day long in the nuclear business takes away a lot
of time that could be used in the "doing without thinking" part. Thus, a
clear definition of the design basis---in the algebraic sense of tiny
vectors that span an arbitrary huge space---was done at very early
stages of the development. Actually, milonga's design basis that follows
was presented at the 2010 Annual Meeting of the Argentine Nuclear
Technology Association [@milonga2010].


There are four main subspaces in the design basis, each spanned by
several vectors, as schematically illustrated in
figure [\[fig:basis\]](#fig:basis){reference-type="ref"
reference="fig:basis"}. The first subspace is about the kind of problems
to solve, the second is about what kind of features the code should have
and should be able to handle, the third is about the expected results
and how to present them to the user and the fourth is about scalability.
Event though they are somehow related to each other, they are discussed
separately in the next four paragraphs. The construction "should be able
to" is deliberately over-used.

![[\[fig:basis\]]{#fig:basis label="fig:basis"}Design basis vectors,
spanning four subspaces](basis)


For the code to be of interest, it should be able to solve detailed
models of both power and research nuclear reactors. So this is one of
the main vectors to keep in mind: the ability to cope with mathematical
descriptions of real cases, incorporating means to take into account the
influence of each of the different parameters that define the actual
flux distribution that real reactors have. Nevertheless, milonga is part
of a design optimization suite and as such, most of the cases in
practice will be conceptual ideas or very crude simplifications of the
final to-be-designed reactor. Thus, simplified cases in one or two
dimensions should also be handled. Besides, because of its academic
nature, the code should also be able to solve cases with analytical
solutions to benchmark numerical schemes and solutions methods. In
addition, the parameters to be optimized usually change the set of
macroscopic cross sections in very diverse ways, so a very flexible way
of providing their dependence on a few to-be-optimized parameters. The
actual approach is discussed in the next subspace, but the important
thing to take into account in this one is that during the optimization,
there may be intermediate steps with parameters that might be
inconsistent or either give rise to unphysical sets of cross sections.
Milonga should be able to cope with these situations without crashing
catastrophically. And last but not least, a great deal of the design
optimization process is based on parametric studies, i.e. analyzing how
certain figures or functions change with a certain parameter while
keeping constant the rest. Solving this parametric kind of problems
should be also a central part of the design basis vectors.


The second subspace spanned by the basis is about flexibility. It is the
main subject around milonga's design and---at least in these firsts
versions---flexibility should have precedence over efficiency. One
important part of this subspace is about the way input data is entered.
A large amount of engineering codes still rely on the card concept, that
is anachronistic, obfuscated, makes no sense nowadays and renders the
preparation of the problem a complex and time-consuming task without
adding extra value. An input preparation concept based in a parser
similar to how compilers translate a human-readable source code into a
binary machine-readable executable should be preferred. Milonga should
read one or more text files containing keywords and arguments that
should be parsed and subsequently converted into the proper coefficients
in the neutron diffusion equation. Another aspect of the flexibility
that is desired in modern engineering computer codes, and in the same
direction as that of the parametric vector, is that the reactor geometry
and the parameters spatial distribution should be defined independently
from the particular spatial nodalization chosen to discretize the
diffusion problem. The spatial distribution of cross sections that
finally characterize the problem to be solved should be viewed as
continuous multidimensional functions with an arbitrary dependence on
other parameters, such as temperature or burn-up distributions that, in
turn, are also viewed as continuous spatial distributions. To handle
varying parameters, milonga should be able to deal with auxiliary
variables, vectors and multidimensional functions, to operate by
applying algebraic or differential operators to them, and to use these
results to evaluate either the geometry definition or cross sections
distributions. The basic results computed by the code should also be
prone to further mathematical manipulation to show only what the user
requests, trying to minimize or avoid the necessity of external data
processing. One big deal of the degree of flexibility of this kind of
computational code is the set of possible sources for the information
that fully define the problem (i.e. geometry, including the position of
the control rods, temperatures, burn-up and poisons distributions). They
might be given as point-wise defined functions or as algebraic
expressions. From step to step, either the independent (the parameters
themselves) or the dependent (their locations) values may change. And of
course, even both. They might be entered in the input, read from local
or remote text files or exchanged with other engineering codes using
some efficient coupling mechanism.


[\[t-rex\]]{#t-rex label="t-rex"} Back when input data was entered by
making holes in a card, T-Rexes ruled the continent and calculation
times were measured in weeks, computer simulation codes had to give as
much output as possible to reduce costs. In the engineering departments,
it was preferred to have a large number of sheets with a lot of tiny
scientific-notation matrix-dot printed numbers stored in cabinets and
shelves for eventual consultation, than to have to re-run simulations
each time a particular result was needed. This is no longer true,
especially for engineering design calculations. Nowadays, the activity
of browsing through old-fashioned huge text files looking for a needle
or having to convert and process numerical data to feed graphical
plotters or post-processing tools is usually far more time consuming
that the execution of the actual simulations. Moreover, as computational
capacity has grown exponentially over the last years---and of course it
is expected to continue increasing---more and more detailed models are
being utilized, and therefore more and more information can be computed.
If all the results are written to the output---consider transient cases
of several days of operation for example---the simulation turns into an
inefficient process in terms of data storage and retrieval. Milonga
should be as flexible as possible in terms of what its output is. First,
no unwanted information should be obtained. Second, output
post-processing and further treatment is to be reduced or even avoided.
And third, means to easily compare results---either with the same
calculation method but with different parameters or with the same
parameters but with different methods---should be taken into account in
the design. Not only should the actual output routines comply with the
basis discussed in this paragraph, but also the whole code structure has
to provide the needed flexibility to present to the user the proper
output information.


Finally, a very important subset of the design basis is that of
scalability. This concept basically means that the code should be able
to run efficiently not only in what are considered current scientific
computing standard architectures, but also to be smoothly adaptable and
take advantage of future improvements in hardware development without
the necessity to be rewritten from scratch. In particular, the
efficiency of an engineering code is directly related to the numerical
methods routines. In milonga, any relevant numerical calculation should
be implemented by using existing software libraries as, for sure,
mathematicians and computer science professionals write better numerical
methods routines than I do. Moreover, relying on well-known
state-of-the-art numerical libraries guarantees scalability, at least up
to the same level the library scales. Also, scalability is closed
related to portability, as this feature increases the chances that the
code could be compiled in whatever architecture is going to be
considered standard for scientific computing in a reasonable software
lifetime (a decade or so). This can only be achieved on the one hand by
sticking to reliable programming standards as much as possible and, on
the other hand, delegating non-specific tasks to available libraries
that also incorporate portability and scalability in their design basis.
In addition, being an academic project, it would be desirable to have a
platform where to test new numerical methods or to compare performances.
Therefore, a straightforward way of incorporating user-coded numerical
methods should be provided. Usually, engineering codes claim to be
modular in nature, easily allowing the incorporation of new features.
But, more often than not, the truth is that they are not modular in the
sense of dynamically-loaded modules that extend the functionality of a
certain computer program---as for example Linux kernel modules do---but
the modularity comes from the use of so-called Fortran modules. This
particularity does not provide the expected benefits of a modular
design---as the Linux kernel does---and, moreover, the overall coding
scheme is prone to obfuscation. Also, the same functionality can be
implemented in a cleaner way using other data structures and languages
like C instead of Fortran, so modularity is not part of milonga's design
basis. To help milonga to survive to changes in computing paradigms and
to scale its power along future hardware availability, its source code
should be freely available for anyone to be able to modify it as
desired, including changes due to flexibility, efficiency, portability
and scalability. Of course, being milonga free software, it can only
depend on free libraries, that may themselves be modified to enhance
flexibility, efficiency, portability and scalability, closing a positive
feedback loop where everyone benefits from the freedom that free
software provides.

Software license {#sec:lincense}
----------------

Living in a country that should be ashamed---amongst other
things---about its software piracy rate [@bsa2010], it seems appropriate
at this point to clearly state the license under which the code is
distributed. To know what kind of rights you are entitled to is as
important as to know how to use the software.


Milonga is free software---both as in free speech and as in free beer,
although the first meaning is far more important than the second
one---and is distributed under the terms of the GNU General Public
License version 3 [@gpl]. In words of the Free Software Foundation,

> Nobody should be restricted by the software they use. There are four
> freedoms that every user should have:
>
> 1.  the freedom to use the software for any purpose,
>
> 2.  the freedom to change the software to suit your needs,
>
> 3.  the freedom to share the software with your friends and neighbors,
>     and
>
> 4.  the freedom to share the changes you make.
>
> When a program offers users all of these freedoms, we call it free
> software.
>
> Developers who write software can release it under the terms of the
> GNU GPL. When they do, it will be free software and stay free
> software, no matter who changes or distributes the program. We call
> this copyleft: the software is copyrighted, but instead of using those
> rights to restrict users like proprietary software does, we use them
> to ensure that every user has freedom.

Not only does milonga provide all the four basic freedoms to the
software user, but also encourages her to study, understand, analyze and
hack it. And of course, to share under the terms of the GNU
GPL---especially with milonga's original author---her discoveries,
suggestions, improvements and fixed bugs. To sum up:


> Milonga is free software: you can redistribute it and/or modify\
> it under the terms of the GNU General Public License as published by\
> the Free Software Foundation, either version 3 of the License, or\
> (at your option) any later version.
>
> 
> Milonga is distributed in the hope that it will be useful,\
> but WITHOUT ANY WARRANTY; without even the implied warranty of\
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\
> GNU General Public License for more details.
>
> 
> You should have received a copy of the GNU General Public License\
> along with wasora. If not, see <http://www.gnu.org/licenses/>.

Milonga relies on a few libraries, all of them available also under
different free licenses:

-   PETSc (<http://www.mcs.anl.gov/petsc>), released freely under the
    copyright of University of Chicago

-   SLEPc (<http://www.grycap.upv.es/slepc>), released freely under the
    GNU Lesser General Public License v3

-   GSL (<http://www.gnu.org/software/gsl>), released freely under the
    GNU General Public License v3

-   Cubature (<http://ab-initio.mit.edu/wiki/index.php/Cubature>),
    released under the GNU General Public License v2 or later

These libraries may depend on other free libraries themselves. Check the
associated documentation for more information.

\bibliographystyle{unsrt}
The equations inside milonga {#cap:equations}
============================


\hspace{\fill}
\sf 
\footnotesize
I am not scared by this complicated scattering\
kernel integral, but by the tiny divergence term.


*Javier Fernandez, PhD in Mathematics,\
after his first looking at the diffusion equation, 2006*


Milonga is a computer code that essentially solves a certain
mathematical equation whose solution, in some way, can be useful from an
engineering point of view. Needless to say, the code is by no means
immune to the garbage in-garbage out concept. Therefore, for the results
to be of real interest, the user should be aware of the mathematics that
milonga is based on, and ought not to execute the code as a black box.
Moreover, the input allows certain options that can be only understood
by knowing what equations milonga solves and what models and
simplifications are implied within them. And last but not least, the
uncertainties and errors introduced by the discretization of a
continuous equation into the finite number of unknowns depend on how the
input continuous parameters are also condensed into a finite set. All
these features are thoroughly explained in this chapter and thus, it is
of central importance for the practical usage of the software.


Milonga is a neutron diffusion code, i.e. it solves the
steady-state---at least in this version---neutron diffusion equation.
This means that its results cannot be more accurate than the diffusion
equation itself, that is already an approximation to the neutron
transport problem---which might be thought of as another simplification
of the physical problem also. Indeed,
figure [\[fig:uncertainties\]](#fig:uncertainties){reference-type="ref"
reference="fig:uncertainties"} shows one of the many conceptual paths
that could be taken to go from the actual real situation to the results
obtained by using a digital computer to solve the problem. Each
rectangle is a source of uncertainties that should be taken into
account, and whose magnitude the user ought to be able to quantify when
executing the last step shown in the diagram.


\reversemarginpar
Whenever an approximation is introduced in the development, an
exclamation mark as the one shown in the left margin will appear.


This chapter focuses on the particular mathematical steps that milonga
takes to go from the continuous diffusion equation up to the numerical
solution of the eigenvalue problem, specifically the three rectangles
shaded in
figure [\[fig:uncertainties\]](#fig:uncertainties){reference-type="ref"
reference="fig:uncertainties"}. Naturally, corrections, suggestions and
improvements are more than welcome. Indeed, the usual \$2.56 reward for
each bug found [@knuth] can also be discussed.

![[\[fig:uncertainties\]]{#fig:uncertainties
label="fig:uncertainties"}One of the many possible paths that connect
the physical problem with the results obtained in a computer. This
chapter focuses just in the route that goes from the diffusion equation
up to---and including---its numerical solution.](uncertainties.pdf)

The neutron diffusion equation
------------------------------

The derivation of the diffusion equation from the neutron transport
theory can be found in the nuclear reactor theory literature. Classical
books that treat the subject paying special attention to mathematical
steps include [@beckurts] and [@glasstone], while other are based on a
physical background like [@henry] and [@duderstadt]. A modern approach
is given in [@cacuci-prinja]. Reference [@lamarsh section 5-3 page 129]
gives explicitly the conditions where the diffusion approximation holds
which, of course, should be met for milonga to be useful. It is assumed
that the user thoroughly understands the physics behind Fick's law and
what the diffusion equation implies.

 
The main goal is to solve the associated critical reactor neutron
diffusion equation with fission sources, namely

![[\[fig:domain\]]{#fig:domain label="fig:domain"}Four-dimensional phase
space where the solution of the diffusion equation is sought for. Three
coordinates correspond to a spatial domain and the other one corresponds
to neutron energy.](domain){width="12cm"}

$$\begin{aligned}
\label{eq:diffusion-continuous-P}
 0 &= \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, E, \mathcal{P}) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},E)\big]\; \Big]\;
- \Sigma_t(\ensuremath\mathbf{r},E,\mathcal{P}) \cdot \phi(\ensuremath\mathbf{r},E) \nonumber \\
 & \quad\quad + \, \int_0^{\infty} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E, \mathcal{P}) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime \nonumber \\
 & \quad\quad + \, \chi(E) \int_0^\infty \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime,\mathcal{P})}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime\end{aligned}$$
for the steady-state flux distribution $\phi(\ensuremath\mathbf{r},E)$
and the effective multiplication factor $k_\text{eff}$ over a
$m+1$-dimensional phase space with spatial coordinates
$\ensuremath\mathbf{r} \in \mathbb{R}^m$ and one energy component $E$,
as depicted in figure [\[fig:domain\]](#fig:domain){reference-type="ref"
reference="fig:domain"}. The divergence and gradient operators are
applied only to the spatial coordinates, and thus boundary
conditions---either Dirichlet, Neumann or mixed---are needed only in the
spatial domain boundary. $D(\ensuremath\mathbf{r}, E, \mathcal{P})$,
$\Sigma_t(\ensuremath\mathbf{r},E,\mathcal{P})$ and
$\nu \Sigma_f(\ensuremath\mathbf{r},E,\mathcal{P})$ are the diffusion
coefficient, total and $\nu$-fission macroscopic cross sections at the
spatial position $\ensuremath\mathbf{r}$ for neutrons of energy $E$.
They can also depend on arbitrary parameters $\mathcal{P}$ that are
further discussed below. The scattering kernel is such
that $\Sigma_s(\ensuremath\mathbf{r},E\rightarrow E^\prime)\,dE\,dE^{\prime}$
gives the macroscopic cross-section in $\ensuremath\mathbf{r}$ for a
scattering collision of a neutron that changes its energy from $E+ dE$
to $E^\prime + dE^\prime$.
Equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"} assumes that there is only one
fissile isotope---or equivalently, that all the fissile isotopes produce
fission neutrons with the same energy distribution, that is indeed the
most common situation. The fission spectrum gives the
probability $\chi(E)\, dE$ for a fission neutron to be born with an
energy between $E$ and $E + dE$. For simplicity, this chapter assumes
there is a single common fission spectrum.


The typical problem to be tackled by milonga is a three-dimensional
reflected nuclear reactor core, including its reactivity control
mechanisms and possibly taking into account the effects of temperature,
density and poison spatial distributions as depicted in
figure [\[fig:typical\]](#fig:typical){reference-type="ref"
reference="fig:typical"}. These distributions can be entered by means of
multidimensional tables, algebraic expressions or read from shared
memory objects (i.e. computed by other codes coupled to milonga). By
using either interpolation or algebraic parsing techniques, the
parameters and thus the nuclear cross sections can be evaluated---and of
course integrated and differentiated---at any point in space as desired.

![[\[fig:typical\]]{#fig:typical label="fig:typical"} Typical problem to
be tackled with milonga. The core and reflector may have different
temperatures, void fractions, poisons, burn-up distribution, etc. These
parameters can be computed by other codes coupled to milonga or they can
be given from the input as multidimensional tables, algebraic
expressions or a combination of all of them. As a result, all the
parameters are continuous functions that can be evaluated at any point
in space.](typical){height="8cm"}


[\[par:nonlinear\]]{#par:nonlinear label="par:nonlinear"} If the
macroscopic nuclear parameters $D$, $\Sigma_t$, $\Sigma_s$ and
$\nu \Sigma_f$ are known functions of the space and the energy only,
then
equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"} is linear. However, in all
practical cases the nuclear properties of the materials depend on other
arbitrary parameters $\mathcal{P}$ that may themselves depend either on
the actual flux at a particular point of the phase space such as fuel,
coolant or moderator temperature, neutronic poison concentration, etc.
In this case,
equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"} is nonlinear.

Moreover, even though the dependence of the macroscopic cross sections
on these parameters can be known *a priori*, the calculation of the
actual value of the parameters is out of milonga's scope and should be
handled by thermal-hydraulic or control system engineering codes.
Nevertheless, this arbitrary dependence can implemented by means of
successive coupled calculations until an certain flux
distribution $\phi^\star(\ensuremath\mathbf{r},E)$ is obtained such that
the corresponding parameters $\mathcal{P}^\star$ give rise to the same
flux $\phi^\star(\ensuremath\mathbf{r},E)$ when inserted back into
equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"}

$$\label{eq:fixed-point}
 \phi^\star(\ensuremath\mathbf{r},E) \rightarrow \mathcal{P}^\star \rightarrow \phi^\star(\ensuremath\mathbf{r},E)$$

For example, consider the case of a reactor calculation involving the
effects of xenon poisoning. The vector of nuclear parameters
$\ensuremath\mathbf{N}$ at position $\ensuremath\mathbf{r}$ depend on
the $^{135}$Xe concentration at $\ensuremath\mathbf{r}$, that depends on
the local neutron flux at $\ensuremath\mathbf{r}$ that of course depends
on the nuclear parameters at $\ensuremath\mathbf{r}$. To solve this
nonlinear problem, first guess an initial xenon
concentration $X_0(\ensuremath\mathbf{r})$---which may be identically
zero. Now evaluate the nuclear parameters
$\ensuremath\mathbf{N}(\ensuremath\mathbf{r},E)$, obtain a flux
distribution $\phi_1(\ensuremath\mathbf{r},E)$ and compute again the
associated equilibrium xenon concentration $X_1(\ensuremath\mathbf{r})$.
Use this distribution to evaluate new nuclear parameters, calculate the
flux and so on.

In some way, a mapping is defined such that at step $n+1$ one has

$$\label{eq:iteration}
 \phi_{n+1} = \ensuremath\mathbf{F}( \phi_{n} )$$ and the actual
solution to the non-linear problem is obtained whenever this mapping has
an stable fixed point that attracts nearby phase-space solutions, as
assumed in
equation [\[eq:fixed-point\]](#eq:fixed-point){reference-type="eqref"
reference="eq:fixed-point"}. The fixed point is obtained when

$$\label{eq:iterations_convergence}
 \phi_{n+1}=\phi_n$$

The fact that there exists a fixed point and that it is stable depends
on $\ensuremath\mathbf{F}(\phi)$, whose behavior is quite difficult to
analyze. However, under normal circumstances---i.e. usual dependence of
cross sections with parameters---there is an stable fixed point that is
the actual solution to the non-linear problem.


In this example case, $\mathcal{P}^\star$ was related to the equilibrium
xenon concentration, but the same idea holds for other parameters such
as temperature and density distributions, boron concentration, etc. Note
that these distributions are outside milonga's scope so in order to take
these effects into account a coupled calculation between the reactor and
other codes is needed.


For a definite step $n$ of this iterative procedure with fixed
parameters, the explicit dependence on $\mathcal{P}$ may be dropped and
thought of as implicit in the dependence of the cross sections on the
position $\ensuremath\mathbf{r}$ and the diffusion equation may be
written as

$$\begin{aligned}
\label{eq:diffusion-continuous}
 0 &= \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, E) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},E)\big]\; \Big]\;
- \Sigma_t(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \nonumber \\
 & \quad\quad + \, \int_0^{\infty} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime 
  + \, \chi(E) \int_0^\infty \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime\end{aligned}$$
that is a linear and homogeneous integro-differential equation.

### Multigroup formulation {#sec:multigroup}

The diffusion
equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"} spans several orders of magnitude
in energy---typically from $10^{-2}$ up to $10^{6}$ eV---with the
macroscopic cross sections varying greatly and abruptly also by
different orders because of nuclear resonances. To handle this, milonga
works with the common multi-group formulation of the diffusion problem
[@henry] which is introduced in this section.

The continuous energy dependence of both the neutron flux and the
nuclear parameters can be transformed into a set of equations describing
the average behavior of the neutrons inside finite intervals of energy
by using the following approach. The continuous energy domain is divided
into $G$ groups with cut-off
values $0=E_G < E_{G-1} < \dots < E_1 < E_0$ not necessarily
equally-spaced, as depicted in
figure [\[fig:multigroup-energy\]](#fig:multigroup-energy){reference-type="ref"
reference="fig:multigroup-energy"}. Energy $E_0$ should be greater than
the energy of the fastest expected neutron. Energy group $g$ is defined
as the interval $[E_{g}, E_{g-1}]$.

![[\[fig:multigroup-energy\]]{#fig:multigroup-energy
label="fig:multigroup-energy"}Partition of the energy range $0<E<E_0$
into $G$ discrete groups](multigroup-energy)

The flux in group $g \in \mathbb{N} \leq G$ is defined as

$$\label{eq:energy-integrated-flux}
 \phi_g(\ensuremath\mathbf{r}) = \phi(\ensuremath\mathbf{r},g) = \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE$$

To avoid an excessive use of sub-indexes, the notation
$\phi(\ensuremath\mathbf{r},g)$ will be preferred. In the analysis that
follows, functions of the integer index $g$ are group values. Note that
the flux in group $g$ is an integrated flux, while the continuous flux
$\phi(\ensuremath\mathbf{r},E)$ is a density in energy. The former has
units of inverse squared length and inverse time (commonly
$\text{cm}^{-2} \, \text{s}^{-1}$) while the latter has units of inverse
squared length, inverse time and inverse energy (i.e.
$\text{cm}^{-2} \, \text{s}^{-1} \, \text{eV}^{-1}$).


The multi-group formulation seeks to obtain $G$ equations for the group
fluxes. Integrating
equation [\[eq:diffusion-continuous\]](#eq:diffusion-continuous){reference-type="eqref"
reference="eq:diffusion-continuous"} with respect to energy between
$E_g$ and $E_{g-1}$ it is

$$\begin{aligned}
 0 &= \int_{E_g}^{E_{g-1}} \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, E) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},E)\big]\; \Big]\; \, dE
- \int_{E_g}^{E_{g-1}} \Sigma_t(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE \nonumber \\
 & \quad\quad + \, \int_{E_g}^{E_{g-1}} \int_0^{\infty} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime  \, dE \\
 & \quad\quad + \, \int_{E_g}^{E_{g-1}} \chi(E) \int_0^\infty \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime \, dE\end{aligned}$$

Now, it is desired to express each term as a product of a group-mean
parameter times a group flux. Starting with the total removal
term---which is the easiest one---one would like to write

$$\int_{E_g}^{E_{g-1}} \Sigma_t(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE = \Sigma_t(\ensuremath\mathbf{r},g) \cdot  \phi(\ensuremath\mathbf{r},g)$$
so the mean total cross section in group $g$ should be defined as

$$\label{eq:sigma_t_weighted}
 \Sigma_t(\ensuremath\mathbf{r},g) = \frac{\displaystyle  \int_{E_g}^{E_{g-1}} \Sigma_t(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE}$$

Of course, if the mean group cross section depends on the flux itself,
then the proposed integration operation is of no help at all. However,
lattice codes are able to obtain reasonable nuclear parameters for the
multi-group formulation by assuming certain energy distributions for the
neutron flux. Indeed, if the flux is assumed to be constant for the
whole energy interval the total cross section is just the traditional
mean value

$$\label{eq:sigma_t_nonweighted}
 \Sigma_t(\ensuremath\mathbf{r},g) = \frac{\displaystyle  \int_{E_g}^{E_{g-1}} \Sigma_t(\ensuremath\mathbf{r},E) \, dE}{\displaystyle E_{g-1} - {E_g}}$$

From milonga's point of view, the value of the macroscopic absorption
cross section is a known function of the spatial coordinates and
eventually of other parameters than can be evaluated before solving the
diffusion equation, at least for each step of the non-linear iteration
[\[eq:iteration\]](#eq:iteration){reference-type="eqref"
reference="eq:iteration"}. Whether the lattice code utilized to generate
the macroscopic cross sections uses
equation [\[eq:sigma\_t\_weighted\]](#eq:sigma_t_weighted){reference-type="eqref"
reference="eq:sigma_t_weighted"} with a certain approximation for
$\phi(\ensuremath\mathbf{r}, E)$,
equation [\[eq:sigma\_t\_nonweighted\]](#eq:sigma_t_nonweighted){reference-type="eqref"
reference="eq:sigma_t_nonweighted"} or even another approach, does not
depend on milonga. The user should know that the code expects a total
macroscopic cross section whose meaning should be compatible with the
definition given by
equation [\[eq:sigma\_t\_weighted\]](#eq:sigma_t_weighted){reference-type="eqref"
reference="eq:sigma_t_weighted"}.


The scattering term can be written as

$$\begin{aligned}
\int_{E_g}^{E_{g-1}} \int_0^{\infty} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime  \, dE &= \\
\int_{E_g}^{E_{g-1}} \sum_{g^\prime = 1}^G \int_{E_{g^\prime}}^{E_{g^\prime-1}} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime  \, dE &= \\
\sum_{g^\prime = 1}^G \int_{E_g}^{E_{g-1}} \int_{E_{g^\prime}}^{E_{g^\prime-1}} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime  \, dE & = \sum_{g' = 1}^G \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g) \\\end{aligned}$$

As the energies $E_g$ are arbitrary, for this equality to hold the
individual terms should be equal. Thus the scattering cross section from
group $g^\prime$ to group $g$ is

$$\label{eq:sigma_s_weighted}
 \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) = \frac{\displaystyle \int_{E_g}^{E_{g-1}} \int_{E_{g^\prime}}^{E_{g^\prime-1}} \Sigma_s(\ensuremath\mathbf{r},E^\prime \rightarrow E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime  \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE}$$

Again, these $G^2$ parameters are treated as known values prior to the
solution of the diffusion equation.


The fission term can be written as

$$\begin{aligned}
\int_{E_g}^{E_{g-1}} \chi(E) \int_0^\infty \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime \, dE &=\\
\int_{E_g}^{E_{g-1}} \chi(E) \sum_{g^\prime = 1}^G \int_{E_{g^\prime}}^{E_{g^\prime-1}} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime \, dE &=\\
\int_{E_g}^{E_{g-1}} \chi(E) \, dE \, \cdot \, \sum_{g^\prime = 1}^G \int_{E_{g^\prime}}^{E_{g^\prime-1}} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},E^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE^\prime & = \chi(g) \cdot \sum_{g^\prime = 1}^G \frac{\nu \Sigma_f(\ensuremath\mathbf{r},g^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r}, g^\prime)\end{aligned}$$
where the fission spectrum of group $g$ is

$$\label{eq:chi_weighted}
 \chi(g) = \int_{E_g}^{E_{g-1}} \chi(E) \, dE$$ and the mean
$\nu$-fission cross section is

$$\label{eq:sigma_f_weighted}
 \nu \Sigma_f(\ensuremath\mathbf{r}, g) = \frac{\displaystyle \int_{E_g}^{E_{g-1}} \nu\Sigma_f(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E^\prime) \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE}$$

Note that the fission spectrum has to be normalized such that

$$\int_{0}^{\infty} \chi(E) \, dE = 1$$ and thus

$$\label{eq:chi_normalized}
 \sum_{g=1}^{G} \chi(g) = 1$$

\vspace{1cm plus 1cm minus 0.25 cm}
Finally, the leakage term---which is the trickiest---has to be written
as

$$\begin{aligned}
~\label{eq:d_exact}
\int_{E_g}^{E_{g-1}} \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, E) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},E)\big]\; \Big]\; \, dE  &= \textsf{div}\, \Big[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},g)\big]\; \Big]\;\end{aligned}$$

In general this equation does not have an exact solution for the
multigroup diffusion coefficient. However, one workaround is to assume
that the flux can be separated into energy and position

$$\begin{aligned}
 \phi(\ensuremath\mathbf{r},E) = \varphi(\ensuremath\mathbf{r}) \cdot \psi(E)\end{aligned}$$
so that

$$\begin{aligned}
\int_{E_g}^{E_{g-1}} \textsf{div}\, \Big[ D(\ensuremath\mathbf{r}, E) \cdot \textsf{grad}\,\big[\varphi(\ensuremath\mathbf{r}) \psi(E)\big]\; \Big]\; \, dE  &= \textsf{div}\, \Big[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},g)\big]\; \Big]\; \\
\textsf{div}\, \left[ \int_{E_g}^{E_{g-1}} D(\ensuremath\mathbf{r}, E) \cdot \psi(E) \cdot \textsf{grad}\,\big[\varphi(\ensuremath\mathbf{r})\big] \, dE\; \right]\; &= \textsf{div}\, \left[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\left[ \int_{E_g}^{E_{g-1}} \varphi(\ensuremath\mathbf{r})\psi(E) \, dE \right]\; \right]\; \\
\textsf{div}\, \left[ \int_{E_g}^{E_{g-1}} D(\ensuremath\mathbf{r}, E) \cdot \psi(E) \, dE \cdot \textsf{grad}\,\big[\varphi(\ensuremath\mathbf{r})\big]\; \right]\; &= \textsf{div}\, \left[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\Big[ \varphi(\ensuremath\mathbf{r})\; \Big]  \cdot \int_{E_g}^{E_{g-1}} \psi(E) \, dE \right]\; \\\end{aligned}$$
and thus the diffusion coefficient for the group $g$ can be computed as

$$D(\ensuremath\mathbf{r},g) = \frac{\displaystyle \int_{E_g}^{E_{g-1}} D(\ensuremath\mathbf{r}, E) \cdot \psi(E) \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \psi(E) \, dE}$$

Another approach may be to calculate the mean transport cross section
for the group $g$ with the same weighting procedure used in the total
cross section
(equation [\[eq:sigma\_t\_weighted\]](#eq:sigma_t_weighted){reference-type="eqref"
reference="eq:sigma_t_weighted"}) and then compute the diffusion
coefficient from its definition

$$D(\ensuremath\mathbf{r},g) = \frac{1}{3 \Sigma_{\text{tr}}(\ensuremath\mathbf{r},g)} = \frac{1}{3} \cdot \frac{\displaystyle \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \Sigma_\text{tr}(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE}$$

Once again, milonga expects the multigroup parameters as known
distributions of space and eventually other parameters such as
temperatures, densities, burn-up and/or poisons, which themselves depend
on the position. So the multigroup parameters, from milonga's point of
view, are essentially known continuous functions of the
position $\ensuremath\mathbf{r}$. It is the user's responsibility to
generate them from a lattice code or whatever other applicable source
consistently with the expected usage shown in this section.


Collecting these results, the continuous diffusion
equation [\[eq:diffusion-continuous\]](#eq:diffusion-continuous){reference-type="eqref"
reference="eq:diffusion-continuous"} can be discretized in energy as $G$
coupled differential equations in the spatial coordinates for the group
neutron fluxes $\phi(\ensuremath\mathbf{r},g)$

$$\begin{aligned}
\label{eq:diffusion_multigroup}
 0 &= \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},g)\big]\; \Big]\;
- \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r},g) \nonumber \\
 & \quad\quad + \, \sum_{g^\prime=1}^{G} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g^\prime) 
  + \, \chi(g) \sum_{g^\prime=1}^{G} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},g^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},g^\prime)\end{aligned}$$

Provided the nuclear parameters are computed from the
definitions [\[eq:sigma\_t\_weighted\]](#eq:sigma_t_weighted){reference-type="eqref"
reference="eq:sigma_t_weighted"},
[\[eq:sigma\_s\_weighted\]](#eq:sigma_s_weighted){reference-type="eqref"
reference="eq:sigma_s_weighted"},
[\[eq:chi\_weighted\]](#eq:chi_weighted){reference-type="eqref"
reference="eq:chi_weighted"},
[\[eq:sigma\_f\_weighted\]](#eq:sigma_f_weighted){reference-type="eqref"
reference="eq:sigma_f_weighted"}
and [\[eq:d\_exact\]](#eq:d_exact){reference-type="eqref"
reference="eq:d_exact"},
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} is exact and the only loss produced
during the energy discretization process is the detailed dependence of
the flux on $E$ inside each group. Any approximation done during the
generation of the multi-group macroscopic parameters will induce
differences between the continuous linear energy diffusion
equation [\[eq:diffusion-continuous\]](#eq:diffusion-continuous){reference-type="eqref"
reference="eq:diffusion-continuous"} and the energy-discretized
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"}.

### Finite volumes spatial discretization {#sec:volumes}

To solve a partial differential equation in a digital computer, a
certain discretization of the spatial coordinates ought to be done. As
one of milonga's design goals is research and development, the code has
several different approaches to the subject, mainly for comparison and
academic reasons. The code is mainly based on a finite volumes scheme as
these kind of methods are best suitable for conservation problems than
finite differences or finite elements are. Nonetheless, some finite
differences schemes recipes are also provided by milonga that are
introduced in section [2.1.3](#sec:differences){reference-type="ref"
reference="sec:differences"}.


The basis of the finite volumes spatial discretization is the division
of the domain into $N$ adjacent cells, as depicted in
figure [\[fig:division-cells\]](#fig:division-cells){reference-type="ref"
reference="fig:division-cells"}. In principle, cells can be arbitrary,
but to fix ideas rectangular and uniform cells are going to be
considered for the moment being, as shown in
figure [\[fig:cells\]](#fig:cells){reference-type="ref"
reference="fig:cells"}. The flux of the group $g \in \mathbb{N} < G$ in
cell $i \in \mathbb{N} < N$ is defined as

$$\label{eq:flux-in-volumes}
 \phi_g^i = \phi(i,g) = \frac{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \, d^m\ensuremath\mathbf{r}}$$
where $\int_{V_i}$ is to be understood as the integral over the volume
of the $i$-th cell and $d^m\ensuremath\mathbf{r}$ represents the volume
differential of the spatial coordinates in $\mathbb{R}^m$. For $m=1,2,3$
(figures [\[fig:cells\]](#fig:cells){reference-type="ref"
reference="fig:cells"}a,
[\[fig:cells\]](#fig:cells){reference-type="ref" reference="fig:cells"}b
and [\[fig:cells\]](#fig:cells){reference-type="ref"
reference="fig:cells"}c) the uniform volume $V_i$ of the $i$-th cell is
equal to $\Delta x$, $\Delta x \Delta y$ and
$\Delta x \Delta y \Delta z$ respectively.

![[\[fig:division-cells\]]{#fig:division-cells
label="fig:division-cells"}Discretization of a two-dimensional spatial
domain into a finite number of cells.](2d-covering.pdf)

\subfloat[1D]{\includegraphics{cell1d}}
\hspace{1cm}
\subfloat[2D]{\includegraphics{cell2d}}
\hspace{1cm}
\subfloat[3D]{\includegraphics{cell3d}}

Again, the notation $\phi(i,g)$ over $\phi_g^i$ will be preferred to
avoid the excessive use of sub-indexes. The flux of group $g$ in
cell $i$ is a scalar number, not a function. The set of the
fluxes $\phi(i,g)$ for $i=1,\dots,N$ and for $g=1,\dots,G$ is the
numerical solution of the diffusion problem sought for. Note
that $\phi(i,g)$ is the mean value of the continuous
flux $\phi(\ensuremath\mathbf{r},g)$ over the $i$-th cell and
thus $\phi(i,g)$ and $\phi(\ensuremath\mathbf{r},g)$ have the same
units---as opposed to what happens
between $\phi(\ensuremath\mathbf{r},E)$
and $\phi(\ensuremath\mathbf{r},g)$ given by
equation [\[eq:energy-integrated-flux\]](#eq:energy-integrated-flux){reference-type="eqref"
reference="eq:energy-integrated-flux"}---namely, inverse squared length
and inverse time.


It should be remarked again that the mean flux $\phi(i,g)$ that will be
the solution to the problem is the mean value of the continuous flux
distribution over the cell, and not the value of the flux evaluated at
the cell center

$$\phi(i,g) = \frac{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \, d^m\ensuremath\mathbf{r}} \neq \phi(\ensuremath\mathbf{r}_i, g)$$
although they are equal in the order of $\Delta x^2$, so the following
approximation may be taken

$$\label{eq:approx-meanflux-centerflux}
 \phi(i,g) \approx \phi(\ensuremath\mathbf{r}_i, g)$$

\marginpar{\hspace{\fill}\includegraphics[height=0.75cm]{approx}}
throughout the mathematical development that follows.

The multigroup diffusion
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} is now integrated throughout the
volume of cell $i$

$$\begin{aligned}
 0 &= \int_{V_i} \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},g)\big]\; \Big]\; \, d^m\ensuremath\mathbf{r}
- \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r},g) \, d^m\ensuremath\mathbf{r}  \\
 & \quad\quad + \, \int_{V_i} \sum_{g^\prime=1}^{G} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r}
  + \, \int_{V_i} \chi(g) \sum_{g^\prime=1}^{G} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},g^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r}\end{aligned}$$
and each term is analyzed separately. First, the term of the total
interactions should be equal to the product of the mean total cross
section, the mean flux and the cell volume

$$\int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r},g) \, d^m\ensuremath\mathbf{r} = \Sigma_t(i,g) \cdot \phi(i,g) \, \int_{V_i} d^m\ensuremath\mathbf{r}$$
so the mean cross section associated to cell $i$ for total neutron
interaction is

$$\label{eq:xs-total-associated}
 \Sigma_t(i,g) = \frac{\displaystyle \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}$$

Now, this parameter is not an input to milonga but should rather be
computed from the continuous total cross section for group $g$ spatial
dependance $\Sigma_s(\ensuremath\mathbf{r},g)$, that is the actual input
to the code. The total cross section may depend
on $\ensuremath\mathbf{r}$ because of materials interfaces and also
because of potential property changes such as temperatures throughout a
single material. As the continuous flux $\phi(\ensuremath\mathbf{r},g)$
is not known, some approximation ought to be done. Milonga provides a
few different methods for computing the integral in
equation [\[eq:xs-total-associated\]](#eq:xs-total-associated){reference-type="eqref"
reference="eq:xs-total-associated"} avoiding direct references to the
flux $\phi(\ensuremath\mathbf{r},g)$. These methods are discussed in
section [2.2](#sec:xs-association){reference-type="ref"
reference="sec:xs-association"}, and give rise to different expressions
for the numerical mean cross sections that should be taken into account
when analyzing results. The rest of this section explicitly develops the
form of the discrete form of the diffusion equation and gives the
equations that should be satisfied by the mean cell cross sections in
order for the discretization to be a faithful representation of the
continuous equation.


The integrated scattering term is

$$\begin{aligned}
\int_{V_i} \sum_{g^\prime=1}^{G} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r} &=\\
\sum_{g^\prime=1}^{G} \int_{V_i} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r} &= \sum_{g^\prime=1}^{G} \Sigma_s(i, g^\prime \rightarrow g) \cdot \phi(i,g^\prime) \, \int_{V_i} d^m\ensuremath\mathbf{r}\end{aligned}$$
and the mean scattering cross section from group $g^\prime$ to group $g$
for cell $i$ is

$$\label{eq:xs-scattering-associated}
\Sigma_s(i, g^\prime \rightarrow g) = \frac{\displaystyle \int_{V_i} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r}, g^\prime) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g^\prime) \, d^m\ensuremath\mathbf{r}}$$

Similarly, the fission term should be

$$\begin{aligned}
 \int_{V_i} \chi(g) \sum_{g^\prime=1}^{G} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},g^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r} &= \\
 \frac{\chi(g)}{k_\text{eff}} \sum_{g^\prime=1}^{G} \int_{V_i} \nu \Sigma_f(\ensuremath\mathbf{r},g^\prime) \cdot \phi(\ensuremath\mathbf{r},g^\prime) \, d^m\ensuremath\mathbf{r} &= \frac{\chi(g)}{k_\text{eff}} \sum_{g^\prime=1}^{G} \nu \Sigma_f(i,g^\prime) \cdot \phi(i,g^\prime) \, \int_{V_i} d^m\ensuremath\mathbf{r}\end{aligned}$$

Therefore

$$\nu \Sigma_f(i,g) = \frac{\displaystyle \int_{V_i} \nu \Sigma_f(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}$$


The only term remaining is the net leakage out of the cell. It is in the
spatial discretization of this term that the fun of the neutron
diffusion problem pops up, so it is analyzed thoroughly for one, two and
three dimensions. An alternative approach based on the finite
differences method is introduced in
section [2.1.3](#sec:differences){reference-type="ref"
reference="sec:differences"}.

#### One dimension

Consider a one-dimensional problem consisting of a slab of length $a$
spanning the interval $[0,a]$ in the $x$ direction surrounded by vacuum.
The interval is divided into $N$ cells of equal length, such that
the $i$-th cell spans the interval $[x_{i-1}, x_i]$, as depicted in
figure [\[fig:1d-discretization\]](#fig:1d-discretization){reference-type="ref"
reference="fig:1d-discretization"}. To simplify the notation, the
neutron flux will be referred to as $\phi(i)$ without any particular
indication of the energy group. The spatial discretization is the same
for all the groups. Integer values of arguments indicate the value of
the mean flux in the cell (crosses) and half-integer arguments indicate
fluxes evaluated at cell boundaries (squares). The first are the values
computed by milonga and the latter are estimated from the former.

![[\[fig:1d-discretization\]]{#fig:1d-discretization
label="fig:1d-discretization"}Spatial discretization of a
one-dimensional domain. The fluxes with integer indexes (crosses) are
the solution sought for, while the fluxes with half-integer arguments
(squares) are estimations from the integer-indexes
fluxes.](1ddiscretization.pdf)


The integrated leakage term in one dimension is

$$\begin{aligned}
\label{eq:leakage-integrated-1d}
 \int_{V_i} \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r})\big]\; \Big]\; \, d\ensuremath\mathbf{r} &= \int_{x_{i-1}}^{x_i} \frac{\partial}{\partial x} \left[ D(x) \cdot \frac{\partial \phi(x)}{\partial x}\right] \, dx \nonumber \\
& = \left[ D(x) \cdot \frac{\partial \phi(x)}{\partial x}\right]^{x_{i}}_{x_{i-1}} \nonumber  \\
& = D(x_i) \cdot \left. \frac{\partial \phi}{\partial x} \right|_{x=x_i} -
D(x_{i-1}) \cdot \left. \frac{\partial \phi}{\partial x} \right|_{x=x_{i-1}}\end{aligned}$$

It may happen that neither $D(x)$ nor $\partial \phi/\partial x$ are
defined at the cell boundary if a material interface exists
at $x = x_i$. However, the product $D(x)\cdot\partial \phi/\partial x$
is always well-defined, as it represents the net neutron current at $x$.
Moreover, the continuity of the current implies that

$$\label{eq:currents-continuity}
 \lim_{\epsilon \rightarrow 0^+} D(x_i - \epsilon) \cdot \left.\frac{\partial \phi}{\partial x}\right|_{x_i - \epsilon} = \lim_{\epsilon \rightarrow 0^+} D(x_i + \epsilon) \cdot \left.\frac{\partial \phi}{\partial x}\right|_{x_i + \epsilon}$$

![[\[fig:1d-intermediateflux\]]{#fig:1d-intermediateflux
label="fig:1d-intermediateflux"}Estimation of the neutron
flux $\phi(i+1/2)$ at the cell boundary (square) from the mean cell
fluxes $\phi(i)$ and $\phi(i+1)$. (crosses). There may be a material
discontinuity at the boundary $x=x_i$, so $\phi(i+1/2)$ may not be equal
to the average of the cell fluxes.](1dintermediateflux.pdf)

The flux $\phi(i+1/2)$ at the boundary between cells $i$ and $i+1$ is
defined so that a discrete first-order version of
equation [\[eq:currents-continuity\]](#eq:currents-continuity){reference-type="eqref"
reference="eq:currents-continuity"} holds

\marginpar{\hspace{\fill}\includegraphics[height=0.75cm]{approx}}
$$D(x_i-\epsilon) \, \frac{\phi(i+1/2) - \phi(i)}{\frac{1}{2} \Delta x} = D(x_i+\epsilon) \, \frac{\phi(i+1) - \phi(i+1/2)}{\frac{1}{2} \Delta x}$$
for some $0 < \epsilon \ll \Delta x$. The flux $\phi(i+1/2)$ as a
function of the fluxes in cells $i$ and $i+1$ is therefore

$$\label{eq:intermediate-flux}
 \phi(i+1/2) = \frac{D(x_i-\epsilon) \cdot \phi(i) + D(x_i+\epsilon) \cdot \phi(i+1)}{D(x_i-\epsilon) + D(x_i+\epsilon)}$$


Note that if $D(x)$ is continuous, the intermediate fluxes in
equation [\[eq:intermediate-flux\]](#eq:intermediate-flux){reference-type="eqref"
reference="eq:intermediate-flux"} are

$$\phi(i+1/2) = \frac{1}{2} \Big[ \phi(i) + \phi(i+1) \Big]$$

The current $D(x) \cdot \partial \phi/\partial x$ at $x_i$ can thus be
approximated as

\marginpar{\hspace{\fill}\includegraphics[height=0.75cm]{approx}}
$$\begin{aligned}
\label{eq:current-xi}
&  \lim_{\epsilon \rightarrow 0^+} D(x_i - \epsilon) \cdot \frac{\partial \phi(x_i-\epsilon)}{\partial x} \approx \nonumber \\
&\quad\quad D(x_i-\epsilon) \left[ \frac{\phi(i+1/2)-\phi(i)}{\frac{1}{2} \Delta x} \right] = \nonumber \\
&\quad\quad D(x_i-\epsilon) \cdot \frac{2}{\Delta x} \left[ \left( \frac{D(x_i-\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)}-1 \right) \phi(i) + \frac{D(x_i+\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)} \cdot \phi(i+1) \right]\end{aligned}$$

Conversely, the current at $x_i$ from the point of view of cell $i+1$
can be written as

$$\begin{aligned}
\label{eq:current-xi-i+1}
&  \lim_{\epsilon \rightarrow 0^+} D(x_i + \epsilon) \cdot \left.\frac{\partial \phi}{\partial x}\right|_{x_i + \epsilon} \approx \nonumber \\
&\quad\quad D(x_i+\epsilon) \left[ \frac{\phi(i+1)-\phi(i+1/2)}{\frac{1}{2} \Delta x} \right] = \nonumber \\
&\quad\quad D(x_i+\epsilon) \cdot \frac{2}{\Delta x} \left[ \left( 1 - \frac{D(x_i+\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)} \right) \phi(i+1) - \frac{D(x_i-\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)} \cdot \phi(i) \right]\end{aligned}$$

Note again that if $D(x)$ is continuous, the current at $x_i$ is

$$D(x_i) \cdot \frac{\phi(i+1) - \phi(i)}{\Delta x}$$

Equation [\[eq:current-xi-i+1\]](#eq:current-xi-i+1){reference-type="eqref"
reference="eq:current-xi-i+1"} can be used to evaluate the current at
$i-1/2$ from the point of view of cell $i$ by replacing $i$ by $i-1$

$$\begin{aligned}
\label{eq:current-xi-1}
&  \lim_{\epsilon \rightarrow 0^+} D(x_{i-1} + \epsilon) \cdot \left.\frac{\partial \phi}{\partial x}\right|_{x_{i-1} + \epsilon} \approx \nonumber \\
&\quad\quad D(x_{i-1}+\epsilon) \left[ \frac{\phi(i)-\phi(i-1/2)}{\frac{1}{2} \Delta x} \right] = \nonumber \\
&D(x_{i-1}+\epsilon) \frac{2}{\Delta x} \left[ \left( 1 - \frac{D(x_{i-1}+\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)} \right) \phi(i) - \frac{D(x_{i-1}-\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)} \cdot \phi(i-1) \right]\end{aligned}$$

Taking into account
equations [\[eq:current-xi\]](#eq:current-xi){reference-type="eqref"
reference="eq:current-xi"}
and [\[eq:current-xi-1\]](#eq:current-xi-1){reference-type="eqref"
reference="eq:current-xi-1"}, the integrated leakage term of
equation [\[eq:leakage-integrated-1d\]](#eq:leakage-integrated-1d){reference-type="eqref"
reference="eq:leakage-integrated-1d"} can be approximated as

$$\begin{aligned}
& D(x_i) \cdot \left. \frac{\partial \phi}{\partial x} \right|_{x=x_i} -
D(x_{i-1}) \cdot \left. \frac{\partial \phi}{\partial x} \right|_{x=x_{i-1}} \approx &\\
&\quad \frac{2 D(x_i-\epsilon)}{\Delta x} \left[ \left( \frac{D(x_i-\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)}-1 \right) \phi(i) + \frac{D(x_i+\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)} \cdot \phi(i+1) \right] - \\
& \frac{2D(x_{i-1}+\epsilon)}{\Delta x} \left[ \left( 1 - \frac{D(x_{i-1}+\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)} \right) \phi(i) - \frac{D(x_{i-1}-\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)} \cdot \phi(i-1) \right]\end{aligned}$$
which can be cleanly rewritten in terms of the cell fluxes at $i$ and
its two first neighbors

$$\label{eq:1d-integrated-leakage-coefficients}
 \int_{x_{i-1}}^{x_i} \frac{\partial}{\partial x} \left[ D(x) \cdot \frac{\partial \phi(x)}{\partial x}\right] \, dx = C_i^{-} \cdot \phi(i-1) + C_i \cdot \phi(i) + C_i^{+} \cdot \phi(i+1)$$
with

$$\begin{aligned}
\label{eq:cs-1d-volumes-epsilon}
 C_{i} &= - \frac{2}{\Delta x} \left[ \frac{D(x_{i-1}-\epsilon) \cdot D(x_{i-1}+\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)}  + \frac{D(x_i-\epsilon) \cdot D(x_i+\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)} \right] \nonumber  \\
\nonumber  \\
 C_i^{-}  &= \frac{2}{\Delta x} \frac{D(x_{i-1}-\epsilon) \cdot D(x_{i-1}+\epsilon)}{D(x_{i-1}-\epsilon)+D(x_{i-1}+\epsilon)} \nonumber \\
\nonumber \\
 C_i^{+} &= \frac{2}{\Delta x} \frac{D(x_i-\epsilon) \cdot D(x_i+\epsilon)}{D(x_i-\epsilon)+D(x_i+\epsilon)}\end{aligned}$$


This set of coefficients gives good results either for cases with $D(x)$
varying smoothly with $x$ or for discontinuities located exactly at cell
boundaries. When there are material interfaces that do not coincide with
the spatial discretization,
equations [\[eq:cs-1d-volumes-epsilon\]](#eq:cs-1d-volumes-epsilon){reference-type="eqref"
reference="eq:cs-1d-volumes-epsilon"} give rise to poor estimations.
This is because if $D(x)$ is continuous at $x=x_i$, then the
intermediate flux $\phi(i+1/2)$ according to
equation [\[eq:intermediate-flux\]](#eq:intermediate-flux){reference-type="eqref"
reference="eq:intermediate-flux"} is always the average of $\phi(i)$ and
$\phi(i+1)$ and thus the current given by
equation [\[eq:current-xi\]](#eq:current-xi){reference-type="eqref"
reference="eq:current-xi"} does not depend explicitly on the actual
position of the interface, but indirectly via the way of computing the
cell-averaged macroscopic cross sections, which is physically incorrect
or at least inaccurate.


One way of avoiding these unphysical results is to use average values
over each half of the cells for the diffusion coefficients

$$\begin{aligned}
\langle D (x_i-\epsilon)\rangle &= \frac{\displaystyle \int_{x_i-\frac{1}{2} \Delta x}^{x_i} D(x) \, dx }{\frac{1}{2} \Delta x}  \\
\\
\langle D (x_i+\epsilon)\rangle  &= \frac{\displaystyle \int_{x_i}^{x_i+\frac{1}{2} \Delta x} D(x) \, dx }{\frac{1}{2} \Delta x}\end{aligned}$$
instead of the actual values $D(x_i - \epsilon)$ and
$D(x_i + \epsilon)$, as illustrated in
figures [\[fig:minusandmean\]](#fig:minusandmean){reference-type="ref"
reference="fig:minusandmean"}. This way, the current depends
continuously on the material interface position and the solutions---for
example the multiplication factor $k_\text{eff}$---also depend
continuously on this position. Whether $D(x_i \mp \epsilon)$ or
$\langle  D(x_i \mp \epsilon)\rangle$ should be used depends on the
particular problem being solved.

\subfloat[]{\includegraphics[width=7cm]{1ddminusepsilon}}
\hspace{0.5cm}
\subfloat[]{\includegraphics[width=7cm]{1dmean}}

Note finally that if the diffusion coefficient is homogeneous in the
whole domain, i.e. $D(x)=D$, the coefficients $C$ (either using $D$
or $\langle D \rangle$) are

$$\begin{aligned}
 C_{i}    &= - \frac{2D}{\Delta x} \\
 C_i^{-}  &= \frac{D}{\Delta x}\\
 C_i^{+}  &= \frac{D}{\Delta x} \end{aligned}$$ and the integrated
leakage term of
equation [\[eq:1d-integrated-leakage-coefficients\]](#eq:1d-integrated-leakage-coefficients){reference-type="eqref"
reference="eq:1d-integrated-leakage-coefficients"} is reduced to

$$\begin{aligned}
 \int_{x_{i-1}}^{x_i} \frac{\partial}{\partial x} \left[ D \cdot \frac{\partial \phi(x)}{\partial x}\right] \, dx \approx \frac{D}{\Delta x} \Big[ \phi(i-1) - 2 \, \phi(i) + \phi(i+1) \Big]\end{aligned}$$
which is the 3-point stencil for the $D$ times the flux Laplacian
multiplied by the cell volume $\Delta x$.

#### Two dimensions

Consider now a two-dimensional case, where the spatial domain is divided
uniformly into $N$ rectangular cells of length $\Delta x$ and
height $\Delta y$, as in
figure [\[fig:cells\]](#fig:cells){reference-type="ref"
reference="fig:cells"}b. Each cell can be uniquely defined by means of
two integer indexes $i_x$ and $i_y$, such that the $i$-th cell (also
referred to as the $i_x$-$i_y$ cell) spans the geometric place
$x_{i_x-1} < x < x_{i_x}$ and $y_{i_y-1}<y<y_{i_y}$
(figure [\[fig:2d-discretization\]](#fig:2d-discretization){reference-type="ref"
reference="fig:2d-discretization"}). The leakage term for this cell is

![[\[fig:2d-discretization\]]{#fig:2d-discretization
label="fig:2d-discretization"}Two-dimensional spatial discretization.
For each energy group, fluxes are computed at two integer indexes
(crosses) and estimated when one index is half-integer
(squares).](2ddiscretization)

$$\int_{V_{i}} \textsf{div}\,\Big[ D(x,y) \cdot \textsf{grad}\,\phi(x,y)\; \Big]\;\, d^2\ensuremath\mathbf{r}$$
that can be rewritten using the divergence theorem as a surface integral

$$\int_{\partial V_{i}} \Big[ D(x,y) \cdot \textsf{grad}\,\phi(x,y)\; \Big] \cdot \hat{\ensuremath\mathbf{n}}(x,y) \,\, d\left( \partial V_{i} \right)$$
where $\partial V_{i}$ is the $i$-th cell boundary and
$\hat{\ensuremath\mathbf{n}}(x,y)$ is the unitary outward vector normal
to the boundary at $\ensuremath\mathbf{r}=(x,y) \in \partial V_{i}$. For
the cell in
figure [\[fig:2d-discretization\]](#fig:2d-discretization){reference-type="ref"
reference="fig:2d-discretization"}, it is

$$\begin{aligned}
\label{eq:2d-4integrals}
\int_{\partial V_{i}} \Big[ D(\ensuremath\mathbf{r}) \cdot \textsf{grad}\,\phi(\ensuremath\mathbf{r})\; \Big] \cdot \hat{\ensuremath\mathbf{n}}(\ensuremath\mathbf{r}) \,\, d\left( \partial V_{i} \right)
=& + \int_{y_{i_y-1}}^{y_{i_y}} D(x_{i_x}, y)    \cdot \frac{\partial \phi}{\partial x}(x_{i_x}, y) \, dy \nonumber \\
& - \int_{y_{i_y-1}}^{y_{i_y}} D(x_{i_x-1}, y)   \cdot \frac{\partial \phi}{\partial x}(x_{i_x-1}, y) \, dy \nonumber \\
& + \int_{x_{i_x-1}}^{x_i} D(x, y_{i_y})         \cdot \frac{\partial \phi}{\partial y}(x, y_{i_y}) \, dx  \nonumber \\
& - \int_{x_{i_x-1}}^{x_i} D(x, y_{i_y-1})       \cdot \frac{\partial \phi}{\partial y}(x, y_{i_y-1}) \, dx \end{aligned}$$

Even though these four integrals ought to be solved individually, they
all have the same structure and can be reduced to results that somehow
resemble the one-dimensional case. Take for example the first one. This
integral represents the net leakage of neutrons passing outward (i.e. to
the right) through the right boundary of the cell. Because of neutron
conservation, it should be equal to the net incoming (i.e. also to the
right) neutron current through the left boundary in
cell $(i_x+1)$-$i_y$. And, because there may be a material interface at
the boundary, proceeding like in the previous section, it should be

$$\label{eq:2d-current-continuity}
\int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x}-\epsilon, y)    \cdot \frac{\partial \phi}{\partial x}(x_{i_x}-\epsilon, y) \, dy =
\int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x}+\epsilon, y)    \cdot \frac{\partial \phi}{\partial x}(x_{i_x}+\epsilon, y) \, dy$$
for some $0 < \epsilon \ll \Delta x$.


To obtain a finite volumes scheme using only four neighbors, is must be
assumed that the net current can be written as

\marginpar{\hspace{\fill}\includegraphics[height=0.75cm]{approx}}
$$\label{eq:2d-first-integral}
\int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x}-\epsilon, y)    \cdot \frac{\partial \phi}{\partial x}(x_{i_x}-\epsilon, y) \, dy
\approx \frac{\partial \phi}{\partial x} \left(x_{{i_x}}-\frac{\Delta x}{2}-\epsilon, y_{i_y}\right) \cdot \int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x}-\epsilon, y) \, dy$$

This approximation is correct up to order $\Delta y^2$. To simplify the
notation in the development that follows, use the nomenclature defined
in
figure [\[fig:2d-discretization\]](#fig:2d-discretization){reference-type="ref"
reference="fig:2d-discretization"}

$$\begin{aligned}
 \phi(i_x,i_y) &= \phi \left(x_{i_x}-\frac{\Delta x}{2}, y_{i_y}-\frac{\Delta y}{2}\right)\\
 \phi(i_x-1/2,i_y) &= \phi \left(x_{i_x-1}, y_{i_y}-\frac{\Delta y}{2}\right)\\
 \phi(i_x+1/2,i_y) &= \phi \left(x_{i_x},   y_{i_y}-\frac{\Delta y}{2}\right)\\
 \phi(i_x,i_y-1/2) &= \phi \left(x_{i_x}-\frac{\Delta x}{2}, y_{i_y-1}\right)\\
 \phi(i_x,i_y+1/2) &= \phi \left(x_{i_x}-\frac{\Delta x}{2}, y_{i_y-1}\right)\\\end{aligned}$$
as the mean flux in the $i_x$-$i_y$ cell and the mean fluxes in each
side of the rectangle, and let

$$\begin{aligned}
 \hat{D}(x,i_y) &= \int_{y_{{i_y}-1}}^{y_{i_y}} D(x, y) \, dy \\
 \hat{D}(i_x,y) &= \int_{x_{{i_x}-1}}^{x_{i_x}} D(x, y) \, dx \\\end{aligned}$$

Note that $x,y \in \mathbb{R}$ and $i_x,i_y \in \mathbb{N}$, so the
first expression is a function of $x$ for a fixed $i_y$ and the second
is a function of $y$ for a fixed $i_x$.


By using the
approximation [\[eq:2d-first-integral\]](#eq:2d-first-integral){reference-type="eqref"
reference="eq:2d-first-integral"} in
equation [\[eq:2d-current-continuity\]](#eq:2d-current-continuity){reference-type="eqref"
reference="eq:2d-current-continuity"} the current continuity implies the
identity

$$\hat{D}(x_{i_x}-\epsilon, i_y) \cdot \frac{\partial \phi}{\partial y} (x_{i_x}-\epsilon,i_y) = 
 \hat{D}(x_{i_x}+\epsilon, i_y) \cdot \frac{\partial \phi}{\partial y} (x_{i_x}+\epsilon,i_y)$$

In the same spirit of
equation [\[eq:currents-continuity\]](#eq:currents-continuity){reference-type="eqref"
reference="eq:currents-continuity"}, the derivatives are approximated by
half increments such that

\marginpar{\hspace{\fill}\includegraphics[height=0.75cm]{approx}}
$$\hat{D}(x_{i_x}-\epsilon, {i_y}) \cdot \frac{\phi({i_x}+1/2,{i_y}) - \phi({i_x},{i_y})}{\frac{1}{2} \Delta x} = 
 \hat{D}(x_{i_x}+\epsilon, {i_y}) \cdot \frac{\phi({i_x}+1,{i_y}) - \phi({i_x}+1/2,{i_y})}{\frac{1}{2} \Delta x}$$

The intermediate flux $\phi(i_x+1/2,i_y)$ as a function of the cell
fluxes is thus

$$\phi(i_x+1/2,i_y) = \frac{\hat{D}(x_{i_x}-\epsilon,{i_y}) \cdot \phi(i_x,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y) \cdot \phi(i_x, i_y+1)}{\hat{D}(x_{i_x}-\epsilon,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y)}$$
and the flux derivative at the boundary cell $i_x$-$i_y$ is

$$\begin{aligned}
\label{eq:2d-current-xi}
 \frac{\partial \phi}{\partial x}(x_{i_x}-\epsilon,i_y) \approx& \frac{2}{\Delta x} \left[ \left( \frac{\hat{D}(x_{i_x}-\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y)} - 1\right) \cdot \phi(i_x,i_y) \right. \nonumber \\
&\quad\quad \left. + \frac{\hat{D}(x_{i_x-1}+\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y)} \cdot \phi(i_x+1,i_y) \right]\end{aligned}$$

Equation [\[eq:2d-current-xi\]](#eq:2d-current-xi){reference-type="eqref"
reference="eq:2d-current-xi"} should be compared with
equation [\[eq:current-xi\]](#eq:current-xi){reference-type="eqref"
reference="eq:current-xi"}. The form is actually the same, with the
two-dimensional version having expressions integrated in the $y$ axis
for the diffusion coefficients, i.e. hats. The first integral in the
right term of
equation [\[eq:2d-4integrals\]](#eq:2d-4integrals){reference-type="eqref"
reference="eq:2d-4integrals"} can then be written as

$$\begin{aligned}
&\int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x}, y)    \cdot \frac{\partial \phi}{\partial x}(x_{i_x}, y) \, dy
\quad \approx \\
&\quad\quad \frac{2 \hat{D}(x_{i_x}-\epsilon,i_y)}{\Delta x} \left[ \left( \frac{\hat{D}(x_{i_x}-\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y)} - 1\right) \phi(i_x,i_y) \right. \\
& \quad\quad\quad\quad\quad\quad\quad\quad
 \quad\quad \left. + \frac{\hat{D}(x_{i_x}+\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y) + \hat{D}(x_{i_x}+\epsilon,i_y)} \cdot \phi(i_x+1,i_x) \right]\end{aligned}$$

An analog reasoning leads to expressions for the rest of the three
integrals as

$$\begin{aligned}
&\int_{y_{{i_y}-1}}^{y_{i_y}} D(x_{i_x-1}, y) \cdot \frac{\partial \phi}{\partial x}(x_{i_x-1}, y) \, dy 
\quad \approx &\\
&\quad\quad \frac{2 \hat{D}(x_{i_x-1}+\epsilon,{i_y})}{\Delta x} \left[ \left( 1 - \frac{\hat{D}(x_{i_x-1}+\epsilon,{i_y})}{\hat{D}(x_{i_x-1}-\epsilon,{i_y}) + \hat{D}(x_{i_x-1}+\epsilon,{i_y})}\right) \phi(i_x,{i_y}) \right. &\\
& \quad\quad\quad\quad\quad\quad\quad\quad\quad\quad
\left. - \frac{\hat{D}(x_{i_x-1}-\epsilon,{i_y})}{\hat{D}(x_{i_x-1}-\epsilon,{i_y}) + \hat{D}(x_{i_x-1}+\epsilon,{i_y})} \cdot \phi(i_x-1,{i_y}) \right] &\end{aligned}$$

$$\begin{aligned}
&\int_{x_{i_x-1}}^{x_{i_x}} D(x, y_{i_y})     \cdot \frac{\partial \phi}{\partial y}(x, y_{i_y}) \, dx
\quad \approx \\
&\quad\quad \frac{2 \hat{D}(i_x,y_{i_y}-\epsilon)}{\Delta y} \left[ \left( \frac{\hat{D}(i_x,y_{i_y}-\epsilon)}{\hat{D}(i_x,y_{i_y}-\epsilon) + \hat{D}(i_x,y_{i_y}+\epsilon)} - 1\right) \phi(i_x,{i_y}) \right. \\
& \quad\quad\quad\quad\quad\quad\quad\quad\quad\quad
 \left. + \frac{\hat{D}(i_x,y_{i_y}+\epsilon)}{\hat{D}(i_x,y_{i_y}-\epsilon) + \hat{D}(i_x,y_{i_y}+\epsilon)} \cdot \phi(i_x,{i_y}+1) \right]\end{aligned}$$

$$\begin{aligned}
&\int_{x_{i_x-1}}^{x_{i_x}} D(x, y_{{i_y}-1}) \cdot \frac{\partial \phi}{\partial y}(x, y_{{i_y}-1}) \, dx
\quad \approx \\
&\quad\quad \frac{2 \hat{D}(i_x,y_{{i_y}-1}+\epsilon)}{\Delta y} \left[ \left( 1 - \frac{\hat{D}(i_x,y_{{i_y}-1}+\epsilon)}{\hat{D}(i_x,y_{{i_y}-1}-\epsilon) + \hat{D}(i_x,y_{{i_y}-1}+\epsilon)}\right) \phi(i_x,{i_y}) \right. \\
& \quad\quad\quad\quad\quad\quad\quad\quad\quad\quad
 \left. - \frac{\hat{D}(i_x,y_{{i_y}-1}-\epsilon)}{\hat{D}(i_x,y_{{i_y}-1}-\epsilon) + \hat{D}(i_x,y_{{i_y}-1}+\epsilon)} \cdot \phi(i_x,{i_y}-1) \right]\end{aligned}$$

Having derived these expressions, the two-dimensional leakage term in
cell $i_x$-$i_y$ can be written as a linear combination of the cell
fluxes at $i_x$-$i_y$ and its four first neighbors

$$\begin{aligned}
\int_{V_{i}} \textsf{div}\,\Big[ D(x,y) \cdot \textsf{grad}\,\phi(x,y)\; \Big]\;\, d^2\ensuremath\mathbf{r} \,=&
 \, C_{i_x,i_y} \, \phi(i_x,i_y) \, +
C_{i_x,i_y}^{x-} \, \phi(i_x-1,i_y) \, +\,  C_{i_x,i_y}^{x+} \, \phi(i_x+1,i_y) \\
& \quad\quad\quad\quad\quad +  C_{i_x,i_y}^{y-} \, \phi(i_x,i_y-1) \, +\,  C_{i_x,i_y}^{y+} \, \phi(i_x,i_y+1)\end{aligned}$$
with the coefficients

$$\begin{aligned}
\label{eq:cs-2d-volumes-epsilon}
 C_{i_x,i_y} =& - \frac{2}{\Delta x} \left[ \frac{\hat{D}(x_{i_x-1}-\epsilon,i_y) \cdot \hat{D}(x_{i_x-1}+\epsilon,i_y)}{\hat{D}(x_{i_x-1}-\epsilon,i_y)+\hat{D}(x_{i_x-1}+\epsilon,i_y)} +   \frac{\hat{D}(x_{i_x}-\epsilon,i_y) \cdot \hat{D}(x_{i_x}+\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y)+\hat{D}(x_{i_x}+\epsilon,i_y)}  \right] \nonumber  \\
& - \frac{2}{\Delta y} \left[ \frac{\hat{D}({i_x},y_{i_y-1}+\epsilon) \cdot \hat{D}({i_x},y_{i_y-1}+\epsilon)}{\hat{D}({i_x},y_{i_y-1}-\epsilon)+\hat{D}({i_x},y_{i_y-1}+\epsilon)}  +   \frac{\hat{D}({i_x},y_{i_y}-\epsilon) \cdot \hat{D}({i_x},y_{i_y}-\epsilon,i_y)}{\hat{D}({i_x},y_{i_y}-\epsilon)+\hat{D}({i_x},y_{i_y}+\epsilon)} \right] \nonumber
\nonumber \\
 C_{i_x,i_y}^{x-}  =& \frac{2}{\Delta x} \frac{\hat{D}(x_{i_x-1}-\epsilon,i_y) \cdot \hat{D}(x_{i_x-1}+\epsilon,i_y)}{\hat{D}(x_{i_x-1}-\epsilon,i_y)+\hat{D}(x_{i_x-1}+\epsilon,i_y)} \nonumber \\
 C_{i_x,i_y}^{y-}  =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,y_{i_y-1}-\epsilon) \cdot \hat{D}(i_x,y_{i_y-1}+\epsilon)}{\hat{D}(i_x,y_{i_y-1}-\epsilon)+\hat{D}(i_x,y_{i_y-1}+\epsilon)} \nonumber \\
 C_{i_x,i_y}^{x+} =& \frac{2}{\Delta x} \frac{\hat{D}(x_{i_x}-\epsilon,i_y) \cdot \hat{D}(x_{i_x}+\epsilon,i_y)}{\hat{D}(x_{i_x}-\epsilon,i_y)+\hat{D}(x_{i_x}+\epsilon,i_y)} \nonumber \\
 C_{i_x,i_y}^{y+} =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,y_{i_y}-\epsilon) \cdot \hat{D}(i_x,y_{i_y}+\epsilon)}{\hat{D}(i_x,y_{i_y}-\epsilon)+\hat{D}(i_x,y_{i_y}+\epsilon)}\end{aligned}$$

As in the previous section, $D(x,y)$ may have a discontinuity that does
not coincide with a cell border, so instead of using the line integral
along $x_i-\epsilon$ to evaluate the net flux through the right, an
average value using the surface integral over the right half of the cell
can be used. For example, $\hat{D}(x_i - \epsilon, i_y)$ may be replaced
by

$$\left\langle \hat{D}(x_i - \epsilon)\right\rangle = \frac{\displaystyle \int_{x_{i_x}-1/2 \Delta x}^{x_{i_y}} \hat{D}(x, i_y) \, dx }{\frac{1}{2} \Delta x} = \frac{\displaystyle \int_{x_{i_x}-1/2 \Delta x}^{x_{i_y}} \int_{y_{i_y-1}}^{y_{i_y}} D(x,y) \, dy \, dx}{\frac{1}{2} \Delta x}$$
with similar expressions for the rest of the diffusion coefficients in
equations [\[eq:cs-2d-volumes-epsilon\]](#eq:cs-2d-volumes-epsilon){reference-type="eqref"
reference="eq:cs-2d-volumes-epsilon"}.


If the diffusion coefficient is homogeneous over the whole spatial
domain $D(x,y)=D$, then the leakage term in two dimensions is

$$\begin{aligned}
\int_{V_{i}} \textsf{div}\,\Big[ D(x,y) \cdot \textsf{grad}\,\phi(x,y)\; \Big]\;\, d^2\ensuremath\mathbf{r} \quad \approx \quad &
 D \cdot \left[ \Delta y \left( \frac{\phi(i_x-1,i_y) - 2\phi(i_x,i_y) + \phi(i_x+1,i_y)}{\Delta x}\right)  \right.\\
& \quad  \left. + \Delta x \left ( \frac{\phi(i_x,i_y-1) - 2\phi(i_x,i_y) + \phi(i_x,i_y+1)}{\Delta y} \right) \right]\end{aligned}$$
that is the expression of $D \nabla^2 \phi \, \Delta x \Delta y$ written
with the usual 5-point stencil for the Laplacian.

#### Three dimensions

The one-dimensional case was used to introduce the basic idea of the
scheme, while the two dimensions case was used to show the extension to
multiple spatial dimensions. The equations for three-dimensional
problems are presented as a generalization of
equations [\[eq:cs-2d-volumes-epsilon\]](#eq:cs-2d-volumes-epsilon){reference-type="eqref"
reference="eq:cs-2d-volumes-epsilon"}:

$$\begin{aligned}
\int_{V_{i}} \textsf{div}\,\Big[ D(x,y,z) \cdot \textsf{grad}\,\phi(x,y,z)\; \Big]\;\, d^3\ensuremath\mathbf{r} \,=&
 \, C_{i_x,i_y,i_z} \, \phi(i_x,i_y,i_z) \\
& \quad +  C_{i_x,i_y,i_z}^{x-} \, \phi(i_x-1,i_y,i_z) \, +\,  C_{i_x,i_y,i_z}^{x+} \, \phi(i_x+1,i_y,i_z) \\
& \quad +  C_{i_x,i_y,i_z}^{y-} \, \phi(i_x,i_y-1,i_z) \, +\,  C_{i_x,i_y,i_z}^{y+} \, \phi(i_x,i_y+1,i_z) \\
& \quad +  C_{i_x,i_y,i_z}^{z-} \, \phi(i_x,i_y,i_z-1) \, +\,  C_{i_x,i_y,i_z}^{z+} \, \phi(i_x,i_y,i_z+1)\end{aligned}$$

$$\begin{aligned}
\label{eq:cs-3d-volumes-epsilon}
 C_{i_x,i_y} =& - \frac{2}{\Delta x} \left[ \frac{\hat{D}(x_{i_x-1}+\epsilon,i_y,i_z) \cdot \hat{D}(x_{i_x-1}-\epsilon,i_y,i_z)}{\hat{D}(x_{i_x-1}-\epsilon,i_y,i_z)+\hat{D}(x_{i_x-1}+\epsilon,i_y,i_z)} +   \frac{\hat{D}(x_{i_x}-\epsilon,i_y,i_z) \cdot \hat{D}(x_{i_x}+\epsilon,i_y,i_z)}{\hat{D}(x_{i_x}-\epsilon,i_y,i_z)+\hat{D}(x_{i_x}+\epsilon,i_y,i_z)}  \right] \nonumber  \\
& - \frac{2}{\Delta y} \left[ \frac{\hat{D}({i_x},y_{i_y-1}-\epsilon,i_z) \cdot \hat{D}({i_x},y_{i_y-1}+\epsilon,i_z)}{\hat{D}({i_x},y_{i_y-1}-\epsilon,i_z)+\hat{D}({i_x},y_{i_y-1}+\epsilon,i_z)}  +   \frac{\hat{D}({i_x},y_{i_y}-\epsilon,i_z) \cdot \hat{D}({i_x},y_{i_y}-\epsilon,i_z)}{\hat{D}({i_x},y_{i_y}-\epsilon,i_z)+\hat{D}({i_x},y_{i_y}+\epsilon,i_z)} \right] \nonumber \\
& - \frac{2}{\Delta z} \left[ \frac{\hat{D}({i_x},{i_y},z_{i_z-1}-\epsilon) \cdot \hat{D}({i_x},{i_y},z_{i_z-1}+\epsilon)}{\hat{D}({i_x},{i_y},z_{i_{z-1}}-\epsilon)+\hat{D}({i_x},{i_y},z_{i_z-1}+\epsilon)}  +   \frac{\hat{D}({i_x},{i_y},z_{i_z}-\epsilon) \cdot \hat{D}({i_x},{i_y},z_{i_z}+\epsilon)}{\hat{D}({i_x},{i_y},z_{i_z}-\epsilon)+\hat{D}({i_x},{i_y},z_{i_z}+\epsilon)} \right] \nonumber
\nonumber \\
 C_{i_x,i_y,i_z}^{x-}  =& \frac{2}{\Delta x} \frac{\hat{D}(x_{i_x-1}-\epsilon,i_y,i_z) \cdot \hat{D}(x_{i_x-1}+\epsilon,i_y,i_z)}{\hat{D}(x_{i_x-1}-\epsilon,i_y,i_z)+\hat{D}(x_{i_x-1}+\epsilon,i_y,i_z)} \nonumber \\
 C_{i_x,i_y,i_z}^{y-}  =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,y_{i_y-1}-\epsilon,i_z) \cdot \hat{D}(i_x,y_{i_y-1}+\epsilon,i_z)}{\hat{D}(i_x,y_{i_y-1}-\epsilon,i_z)+\hat{D}(i_x,y_{i_y-1}+\epsilon,i_z)} \nonumber \\
 C_{i_x,i_y,i_z}^{z-}  =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,i_y,z_{i_z-1}-\epsilon) \cdot \hat{D}(i_x,i_y,z_{i_z-1}+\epsilon)}{\hat{D}(i_x,i_y,z_{i_z-1}-\epsilon)+\hat{D}(i_x,i_y,z_{i_z-1}+\epsilon)} \nonumber \\
 C_{i_x,i_y,i_z}^{x+} =& \frac{2}{\Delta x} \frac{\hat{D}(x_{i_x}-\epsilon,i_y,i_z) \cdot \hat{D}(x_{i_x}+\epsilon,i_y,i_z)}{\hat{D}(x_{i_x}-\epsilon,i_y,i_z)+\hat{D}(x_{i_x}+\epsilon,i_y,i_z)} \nonumber \\
 C_{i_x,i_y,i_z}^{y+} =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,y_{i_y}-\epsilon,i_z) \cdot \hat{D}(i_x,y_{i_y}+\epsilon,i_z)}{\hat{D}(i_x,y_{i_y}-\epsilon,i_z)+\hat{D}(i_x,y_{i_y}+\epsilon,i_z)} \nonumber \\
 C_{i_x,i_y,i_z}^{z+}  =& \frac{2}{\Delta y} \frac{\hat{D}(i_x,i_y,z_{i_z}-\epsilon) \cdot \hat{D}(i_x,i_y,z_{i_z}+\epsilon)}{\hat{D}(i_x,i_y,z_{i_z}-\epsilon)+\hat{D}(i_x,i_y,z_{i_z}+\epsilon)}\end{aligned}$$

### Finite differences scheme {#sec:differences}

Instead of a finite volumes scheme, an approach based on the finite
differences method can be used to discretize the neutron diffusion
equation. In this case, the spatial domain is also divided in a finite
number $N$ of cells, but instead of working with cell-averaged values,
the flux is computed at the cell centers that are taken as the problem
nodes as depicted in
figure [\[fig:differences\]](#fig:differences){reference-type="ref"
reference="fig:differences"}.

![[\[fig:differences\]]{#fig:differences
label="fig:differences"}Finite-differences approach to the
discretization of the spatial domain introduced in
figure [\[fig:division-cells\]](#fig:division-cells){reference-type="ref"
reference="fig:division-cells"}. The properties and the flux are
evaluated at nodes (filled circles) instead of representing cell
(squares) mean values. For convenience, the nodes are taken as the
finite-volumes cells centers.](2d-covering-differences)


The finite differences scheme consists in replacing the $G$ continuous
multigroup diffusion
equations [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"}

$$\begin{aligned}
 0 &= \textsf{div}\,\Big[ D(\ensuremath\mathbf{r}, g) \cdot \textsf{grad}\,\big[\phi(\ensuremath\mathbf{r},g)\big]\; \Big]\;
- \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r},g) \nonumber \\
 & \quad\quad + \, \sum_{g^\prime=1}^{G} \Sigma_s(\ensuremath\mathbf{r},g^\prime \rightarrow g) \cdot \phi(\ensuremath\mathbf{r},g^\prime) 
  + \, \chi(g) \sum_{g^\prime=1}^{G} \frac{\nu \Sigma_f(\ensuremath\mathbf{r},g^\prime)}{k_\text{eff}} \cdot \phi(\ensuremath\mathbf{r},g^\prime)\end{aligned}$$
by $NG$ equations with the cross sections evaluated in each cell and the
divergence term replaced by a differences approximation. Denote by the
integer argument $i$ the properties evaluated at the center of
the $i$-th cell in an $m$-dimensional spatial domain. As
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} holds for every point in the
spatial domain, then in particular it must hold for the $N$ cell
centers $\ensuremath\mathbf{r}_i$

$$\begin{aligned}
 0 &= \textsf{div}\,\Big[ D(i, g) \cdot \textsf{grad}\,\big[\phi(i,g)\big]\; \Big]\;
- \Sigma_t(i,g) \cdot \phi(i,g) \nonumber \\
 & \quad\quad + \, \sum_{g^\prime=1}^{G} \Sigma_s(i,g^\prime \rightarrow g) \cdot \phi(i,g^\prime) 
  + \, \chi(g) \sum_{g^\prime=1}^{G} \frac{\nu \Sigma_f(i,g^\prime)}{k_\text{eff}} \cdot \phi(i,g^\prime)\end{aligned}$$

The first term can be expanded as

$$\label{eq:diff-grad}
 \textsf{div}\,\Big[ D(i, g) \cdot \textsf{grad}\,\big[\phi(i,g)\big]\; \Big]\; =
D(i,g) \cdot \nabla^2 \phi(i,g) + \textsf{grad}\, \big[ D(i,g)\big] \; \cdot \textsf{grad}\, \big[\phi(i,g)\big] \;$$

The finite differences scheme proposed consists of approximating both
the flux Laplacian and gradient by first-neighbors differences
expressions. In one dimension, it is

$$\begin{aligned}
 \textsf{div}\,\Big[ D(i, g) \, \textsf{grad}\,\big[\phi(i,g)\big]\; \Big]\; \approx&
\quad D(i,g) \cdot \frac{\phi(i+1,g) -2\phi(i,g)+ \phi(i-1,g)}{\Delta x^2}\\
& \quad + \frac{\partial D(i,g)}{\partial x}  \cdot \frac{\phi(i+1,g) - \phi(i-1,g)}{2 \Delta x}\end{aligned}$$

In two dimensions, using the $i_x$-$i_y$ nomenclature introduced in
section [2.1.2](#sec:volumes){reference-type="ref"
reference="sec:volumes"}, the scheme is

$$\begin{aligned}
 \textsf{div}\,\Big[ D(i_x, i_y, g) \cdot \textsf{grad}\,\big[\phi(i_x,i_y,g)\big]\; \Big]\; \approx&
  \quad D(i_x,i_y,g) \, \frac{\phi(i_x+1,i_y,g) -2\phi(i_x,i_y,g)+ \phi(i_x-1,i_y,g)}{\Delta x^2}\\
& \quad + D(i_x,i_y,g) \, \frac{\phi(i_x,i_y+1,g) -2\phi(i_x,i_y,g)+ \phi(i_x,i_y-1,g)}{\Delta y^2}\\
& \quad + \frac{\partial D(i_x,i_y,g)}{\partial x}  \cdot \frac{\phi(i_x+1,i_y,g) - \phi(i_x-1,i_y,g)}{2 \Delta x}\\
& \quad + \frac{\partial D(i_x,i_y,g)}{\partial y}  \cdot \frac{\phi(i_x,i_y+1,g) - \phi(i_x,i_y-1,g)}{2 \Delta y}\end{aligned}$$
and for three dimensions

$$\begin{aligned}
& \textsf{div}\,\Big[ D(i_x, i_y, i_z, g) \cdot \textsf{grad}\,\big[\phi(i_x,i_y,i_z,g)\big]\; \Big]\; \approx &
\\
& \quad\quad D(i_x,i_y,i_z,g) \, \frac{\phi(i_x+1,i_y,i_z,g) -2\phi(i_x,i_y,i_z,g)+ \phi(i_x-1,i_y,i_z,g)}{\Delta x^2} &\\
& \quad\quad + D(i_x,i_y,i_z,g) \, \frac{\phi(i_x,i_y+1,i_z,g) -2\phi(i_x,i_y,i_z,g) + \phi(i_x,i_y-1,i_z,g)}{\Delta y^2} &\\
& \quad\quad + D(i_x,i_y,i_z,g) \, \frac{\phi(i_x,i_y,i_z+1,g) -2\phi(i_x,i_y,i_z,g) + \phi(i_x,i_y,i_z-1,g)}{\Delta z^2} &\\
&\quad\quad + \frac{\partial D(i_x,i_y,i_z,g)}{\partial x}  \cdot \frac{\phi(i_x+1,i_y,i_z,g) - \phi(i_x-1,i_y,i_z,g)}{2 \Delta x}&\\
&\quad\quad + \frac{\partial D(i_x,i_y,i_z,g)}{\partial y}  \cdot \frac{\phi(i_x,i_y+1,i_z,g) - \phi(i_x,i_y-1,i_z,g)}{2 \Delta y}&\\
&\quad\quad + \frac{\partial D(i_x,i_y,i_z,g)}{\partial y}  \cdot \frac{\phi(i_x,i_y,i_z+1,g) - \phi(i_x,i_y,i_z-1,g)}{2 \Delta z}&\end{aligned}$$

Because of the differential nature of the method, only local effects are
retained in this formulation and variations of parameters within cells
are usually damped an even neglected. To somehow take into account these
effects, two things can be done. The first one is to use cell-averaged
properties instead of the values evaluated at the centers (see
section [2.2](#sec:xs-association){reference-type="ref"
reference="sec:xs-association"}), and the second one is to compute the
gradient of the diffusion coefficient also as a finite-difference
approximation instead of using the actual derivative of the continuous
property. All these four combinations can be handled by milonga.


Last but not least, there is one important thing to add. If the
diffusion coefficient $D(\ensuremath\mathbf{r})$ is not differentiable
at some point $\ensuremath\mathbf{r}_d$, then
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} does not hold
at $\ensuremath\mathbf{r}_d$. If this point happens to coincide with a
cell center $\ensuremath\mathbf{r}_i$, the numerical method is not
consistent. Moreover, even if $\ensuremath\mathbf{r}_d$ does not
coincide with a cell center, the fluxes in the neighboring cells will
have unphysical results because the original
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} is not valid. If one wants to find
the solution to the differential equation all in all, current continuity
conditions have to be solved at $\ensuremath\mathbf{r}_d$ instead of the
diffusion equation. The finite volumes scheme proposed in
section [2.1.2](#sec:volumes){reference-type="ref"
reference="sec:volumes"} includes both the diffusion equation and
current continuity in its formulation, so it is expected to give better
results than the finite differences approximation---at least without
using ad-hoc correction factors. For example, if there is a
discontinuity in the distribution of the diffusion coefficient then the
evaluation of $\textsf{grad}\,\phi\;$ may depend on the size of the
cells.Also note that while the finite differences scheme solves the
differential diffusion equation, the finite volumes scheme solves the
diffusion equation integrated over the cell volume.

Cell cross-section computation {#sec:xs-association}
------------------------------

The discretized multigroup equations give a relationship between the
expected neutron flux in each cell for each energy group (the $\phi$'s)
and some parameters that depend on the nuclear properties of the
materials present in the diffusive media (the $\Sigma$'s). Strictly
speaking, these macroscopic cross sections depend on the local flux
distribution as shown by the exact
equation [\[eq:xs-total-associated\]](#eq:xs-total-associated){reference-type="eqref"
reference="eq:xs-total-associated"} derived for the total macroscopic
cross section:

$$\tag{\ref{eq:xs-total-associated}}
 \Sigma_t(i,g) = \frac{\displaystyle \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}$$

However, it is clear that for the formulation to be useful, the
coefficients ought not to depend on the flux itself.


Up to a first order approximation, two approaches can be taken to take
out the fluxes of the cross-section expression. The first is to assume a
flat profile inside each cell. Call
it $\phi(\ensuremath\mathbf{r},g)=\varphi$ for
$\ensuremath\mathbf{r} \in$ cell $i$. Then

$$\label{eq:cell-xs-mean}
 \Sigma_t(i,g) = \frac{\displaystyle \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \varphi \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \varphi \, d^m\ensuremath\mathbf{r}} = \frac{\displaystyle \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \, d^m\ensuremath\mathbf{r}}$$
and the cell total macroscopic cross section is the mean value of the
continuous macroscopic cross section distribution inside the cell.


The second---and rather extreme--approach is to assume that somehow the
flux is concentrated in the center of the
cell $\phi(\ensuremath\mathbf{r},g)=\varphi \cdot \delta(\ensuremath\mathbf{r}_i)$
where $\ensuremath\mathbf{r}_i$ is the location of the $i$-th cell
center. Thus

$$\label{eq:cell-xs-center}
 \Sigma_t(i,g) = \frac{\displaystyle \int_{V_i} \Sigma_t(\ensuremath\mathbf{r},g) \cdot \varphi \, \delta(\ensuremath\mathbf{r_i}) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_{V_i} \varphi \, \delta(\ensuremath\mathbf{r}_i) \, d^m\ensuremath\mathbf{r}} = \Sigma_t(\ensuremath\mathbf{r}_i,g)$$
and the cell cross section is the value taken by the distribution at the
cell center.

Either of these schemes can be chosen to be used by milonga. They give
similar results for cross sections that depend weakly
on $\ensuremath\mathbf{r}$, so
equation [\[eq:cell-xs-center\]](#eq:cell-xs-center){reference-type="eqref"
reference="eq:cell-xs-center"} should be preferred because it requires
less computational effort. However, for discontinuous
properties---especially if they do not coincide with cell
borders---equation [\[eq:cell-xs-mean\]](#eq:cell-xs-mean){reference-type="eqref"
reference="eq:cell-xs-mean"} should be used. If the cross section
distribution is homogeneous inside the cell, both formulations coincide.

Boundary conditions {#sec:boundary}
-------------------

To fully define the diffusion problem, boundary conditions on the
spatial domain must be given. These can be Dirichlet, Neumann or mixed,
i.e. specify the flux, the derivative in the normal direction or a
linear combination of them at the domain boundary. In the discretized
diffusion problem, as illustrated for a two-dimensional case in
figure [\[fig:boundary\]](#fig:boundary){reference-type="ref"
reference="fig:boundary"}, in those cells that are completely contained
in the domain the diffusion equation is solved. In those that are
completely outside the domain, the flux is forced to be null. And in
those that contain the external boundary of the geometry that defines
the neutronic problem, the boundary conditions equations are imposed.

![[\[fig:boundary\]]{#fig:boundary label="fig:boundary"}For the spatial
domain of
figure [\[fig:division-cells\]](#fig:division-cells){reference-type="ref"
reference="fig:division-cells"}, cells that are completely inside the
domain (green) solve the diffusion problem. Cells that contain the
border (orange) impose the problem boundary conditions. The rest of the
cells (black) are forced to have a null flux.](boundary)


In nuclear reactor diffusion theory [@duderstadt], it is usual to make
the flux vanish at a certain distance from the external surfaces called
extrapolation length. However, for common reactor sizes, this length is
negligible enough to assume it is zero. Thus, flux equal to zero at the
external surfaces is the common boundary condition. For symmetric
problems, it is also customary to solve only a part of the problem and
impose mirror boundary conditions on the flux. The boundary conditions
supported by milonga are introduced in the following subsections,
stating the actual discrete equations that implement them numerically.

### External planar surfaces

External planar surfaces include one-dimensional slab boundaries,
two-dimensional straight lines and three-dimensional planar surfaces
where the flux is ought to be equal to zero. In each case, there is a
well-defined outward direction. Without loss of generality, assume it is
the negative $x$ direction, as if a null boundary condition is applied
at $x=0$ in a one-dimensional slab.

Let $i_0$ be the index of the cell that should impose the boundary
condition equation. In general, the external boundary will be at an
algebraic distance $\delta x$ (i.e. $\delta x > 0$ if the surface is at
the right and conversely) from the cell center $x_{i_0}$. The condition
$\phi(x_{i_0}+\delta x) = 0$ is implemented by requiring that a linear
interpolation (or extrapolation for $\delta x < 0$) of the cell fluxes
at $i_0$ and its first neighbor in the inward direction vanishes at
$x=x_{i_0}+\delta x$
(figure [\[fig:bc-planar\]](#fig:bc-planar){reference-type="ref"
reference="fig:bc-planar"}).

![[\[fig:bc-planar\]]{#fig:bc-planar label="fig:bc-planar"} Planar
boundary conditions in the first cell of the $x$ direction. Note that in
this case, the external surface is located between the bare external
surface and the first cell center, so $\delta x <0$. In the case shown,
$i_0=1$](1dplanarbc)

Making use of
approximation [\[eq:approx-meanflux-centerflux\]](#eq:approx-meanflux-centerflux){reference-type="eqref"
reference="eq:approx-meanflux-centerflux"}, further assume that inside
cell $i_0$ the flux can be written as

$$\begin{aligned}
 \phi(x_{i_0} + \delta x) &= \phi(x_{i_0}) + \left. \frac{\partial \phi}{\partial x} \right|_{x_{i_0}} \cdot \delta x \\
 \phi(x_{i_0} + \delta x) &= \phi(i_0) + \frac{\phi(i_0+1)-\phi(i_0)}{\Delta x} \cdot \delta x\end{aligned}$$

The condition $\phi(x_{i_0} + \delta x)=0$ is attained if

$$\phi(i_0) + \left( \frac{\displaystyle \frac{\delta x}{\Delta x}}{1 - \displaystyle  \frac{\delta x}{\Delta x}} \right) \cdot \phi(i_0+1) = 0$$

In particular, if the boundary condition has to be applied at $x=0$,
then $\delta x=-1/2 \Delta x$ and the boundary condition yields

$$\phi(1) = \frac{1}{3} \, \phi(2)$$ meaning that the flux at the first
cell has to be one-third of the flux in the first neighbor in order for
the linear extrapolation to vanish at the boundary of the cell (try to
make a mental picture and see the similar triangles with ratio three to
one).


[\[negative-fluxes\]]{#negative-fluxes label="negative-fluxes"}It is
important to note that if the boundary lies to the right of the cell
center ($\delta x > 0$), then the flux $\phi(i_0)$ has to be negative to
make the linear interpolation between $\phi(i_0)$ and $\phi(i_0+1)$ pass
through zero as shown in
figure [\[fig:bc-planarneg\]](#fig:bc-planarneg){reference-type="ref"
reference="fig:bc-planarneg"}. If this fact is understood as a
mathematical artifact to handle discrete boundary conditions and the
results are properly interpreted, the cell flux can be left negative as
computed. However, fission power is usually computed from the neutron
fluxes and negative powers cannot be allowed based on physical grounds.
Milonga can be told to force negative fluxes to zero so no artificial
power sinks appear, at the expense of slightly changing the form of the
neutron flux near the borders.

![[\[fig:bc-planarneg\]]{#fig:bc-planarneg label="fig:bc-planarneg"} If
the external surface is located between two adjacent nodes, to make the
linear interpolation of the flux pass through zero at the external
surface, one of the fluxes has to be negative. This is a mathematical
artifact to handle discrete boundary conditions, and should be
understood as such.](1dplanarbcneg)


Had the outward direction been the positive $x$ direction, then the
equation to impose the boundary condition would have been

$$\phi(i_0) - \left( \frac{\displaystyle \frac{\delta x}{\Delta x}}{1 + \displaystyle  \frac{\delta x}{\Delta x}} \right) \cdot \phi(i_0-1) = 0$$
with $\delta x > 0$ still meaning that the external surface is located
to the right of cell $i_0$, and conversely. For two and three
dimensions, the idea is extended to all the possible directions.

### Cylindrical surfaces

For cylindrical external surfaces, the only available boundary condition
is null flux. Assume that the cylinder axis is parallel to the $z$ axis,
then

$$\phi(x, y) = 0$$ for $(x,y)$ in the external boundary.

The discrete boundary condition requires a linear expansion of the flux
to be equal to zero at the point in the circumference that lies in the
circle radius that passes through the cell center. Let $(x_0,y_0)$ be
the coordinates of the center of the cell $(i_0, j_0)$ through which the
circumference with center $(x_c,y_c)$ and radius $R$ passes by
(figure [\[fig:2dcirclebc\]](#fig:2dcirclebc){reference-type="ref"
reference="fig:2dcirclebc"}). The point $(x^\star, y^\star)$ at which
the expansion has to vanish is

![[\[fig:2dcirclebc\]]{#fig:2dcirclebc label="fig:2dcirclebc"}Discrete
boundary condition in a cylindrical external surface](2dcirclebc)

$$\begin{aligned}
 x^\star &= x_c \pm R \cos \theta \\
 y^\star &= y_c \pm R \sin \theta \\\end{aligned}$$ with the sign chosen
according to the quadrant $(x_0,y_0)$ is with respect to $(x_c,y_c)$ and

$$\theta = \arctan \left( \frac{|y_0 - y_c|}{|x_0 - y_c|} \right)$$

The linear expansion is a first-order Taylor series around the cell
center

$$\phi(x,y) \approx \phi(x_0, y_0) \, + \, \left.\frac{\partial \phi}{\partial x}\right|_{x_0,y_0} \cdot (x - x_0) \,+\, \left.\frac{\partial \phi}{\partial y}\right|_{x_0,y_0} \cdot (y - y_0)$$

The partial derivatives are approximated as

$$\begin{aligned}
 \left.\frac{\partial \phi}{\partial x}\right|_{x_0,y_0} &\approx \frac{\phi(i_0 \pm 1,j_0) \mp \phi(i_0, j_0)}{\Delta x} \\
 \left.\frac{\partial \phi}{\partial y}\right|_{x_0,y_0} &\approx \frac{\phi(i_0, j_0 \pm 1) \mp \phi(i_0, j_0)}{\Delta y}\end{aligned}$$
where again the signs are selected according to the quadrant with
respect to the circle center the point $(x_0,y_0)$ is. To fix ideas,
assume it is in the first quadrant. Then, the discrete boundary
condition is

$$\begin{aligned}
0 &= \left[ 1 + \frac{x_c-x_0+R\cos\theta}{\Delta x} + \frac{y_c-y_0+R\sin\theta}{\Delta y} \right] \cdot \phi(i_0, j_0) \\
&\quad\quad\quad - \frac{x_c-x_0+R\cos\theta}{\Delta x} \cdot \phi(i_0-1,j_0) - \frac{y_c-y_0+R\sin\theta}{\Delta y} \cdot \phi(i_0, j_0-1)\end{aligned}$$
with similar equations for the other three quadrants.


Note that depending on the relative location of the cell
center $(x_0,y_0)$ and the point $(x^\star, y^\star)$, the boundary
equation may require the flux $\phi(i_0, j_0)$ to be negative as with
the planar surfaces discussed above. Again, this is a mathematical
artifact that hast to be taken into account in the nodalization of a
problem.

### Spherical surfaces

To be done.

### Mirror conditions

Symmetry conditions at a planar surfaces are implemented by requiring a
null partial derivative in the normal direction. Independently from the
actual location of the surface with respect to the cells $i_0$
and $i_0\pm1$, this condition reads

$$\phi(i_0) - \phi(i_0+1) = 0$$ or $$\phi(i_0) - \phi(i_0-1) = 0$$
according to the outward direction of the surface with the mirror
condition.

### Sharp edges

In milonga's context, a sharp edge means a geometric object that belongs
to the external boundary of the spatial domain that has at least two
degrees of dimension less than the domain. For example, in three
dimensions corner points and external surfaces edges are considered
sharp (figure [\[fig:3dsharp\]](#fig:3dsharp){reference-type="ref"
reference="fig:3dsharp"}). In two dimensions, the vertex of the
intersection of two planar external surfaces are sharp edges
(figure [\[fig:2dsharp\]](#fig:2dsharp){reference-type="ref"
reference="fig:2dsharp"}). In one dimension there are no sharp edges.

\subfloat[\label{fig:3dsharp} In three dimensions, sharp edges include corner and surfaces intersections]{\includegraphics{3dsharp}}
\hspace{3cm}
\subfloat[\label{fig:2dsharp} In two dimensions, sharp edges are corners.]{\includegraphics{2dsharp}}
When a cell contains a sharp edge, instead of solving the diffusion
equation a boundary condition equation is applied to it. As there is no
definite outward direction, the condition in milonga is that the flux at
the cell should be equal to the average flux of its non-vacuum
neighbors. This choice of boundary condition works well both with null
and mirror conditions for the planar surfaces.

The eigenvalue problem
----------------------

As assumed throughout this chapter, in the absence of fixed independent
sources---i.e. not associated to fission---the steady-state multigroup
neutron diffusion equation can be written as a generalized eigenvalue
problem. In effect, the $G$ differential
equations [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} have the same dependency on the
multiplicative factor $k_\text{eff}$, and can be casted in matrix form
as

$$\label{eq:multigroup-eigenvalue}
 R^{\star} \cdot \boldsymbol{\varphi} = \frac{1}{k_\text{eff}} F^{\star} \cdot \boldsymbol{\varphi}$$
where $R^\star$ and $F^\star$ are square matrices of size $G \times G$
that contain differential operators and functions of the spatial
coordinate $\ensuremath\mathbf{r}\in\mathbb{R}^m$,
and $\boldsymbol{\varphi}$ is the vector of size $G$

$$\boldsymbol{\varphi} =
\begin{bmatrix}
 \phi(\ensuremath\mathbf{r},1) \\
 \phi(\ensuremath\mathbf{r},2) \\
 \vdots \\
 \phi(\ensuremath\mathbf{r},G) \\
\end{bmatrix}$$

The star notation means that these matrices represent the differential
formulation, while we leave the non-star notation for the
spatial-discretized problem that is milonga's main subject.
Matrices $R^\star$ and $F^\star$ are named after removal and fission,
respectively. The left member of
equation [\[eq:multigroup-eigenvalue\]](#eq:multigroup-eigenvalue){reference-type="eqref"
reference="eq:multigroup-eigenvalue"} gives the net rate of
disappearance of neutrons, while the right member gives the rate of
neutron births by fission. By using the mathematical artifact of
dividing this rate by the effective multiplication factor, an artificial
critical reactor is obtained, with the eigenvector representing the
artificial flux distribution.
Equation [\[eq:multigroup-eigenvalue\]](#eq:multigroup-eigenvalue){reference-type="eqref"
reference="eq:multigroup-eigenvalue"} can be written as a standard
eigenvalue problem as

$$\label{eq:standard-eigenvalue}
 {R^{\star}}^{-1} F^{\star} \cdot \boldsymbol{\varphi} = {k_\text{eff}} \cdot \boldsymbol{\varphi}$$

It is possible to show that $R^{-1}$ exists for any physically real set
of cross sections [@henry], so the eigenvalue problem is well-defined.
Now, in general there are $G$ pairs of eigenvalues and eigenvectors that
satisfy the eigenvalue problem. Even though most of them may be
negative, complex or give rise to negative fluxes, for problems based on
a physical ground it is possible to prove [@henry section 3.3, page 73]
that

1.  There is a unique real positive eigenvalue greater in magnitude than
    any other eigenvalue

2.  All the elements of the eigenvector corresponding to that eigenvalue
    are real and positive

3.  All other eigenvectors either have some elements that are zero or
    have elements that differ in sign from each other


The eigenvalue
form [\[eq:multigroup-eigenvalue\]](#eq:multigroup-eigenvalue){reference-type="eqref"
reference="eq:multigroup-eigenvalue"} is maintained for the
spatial-discretized formulation derived in
sections [\[sec:volumes\]](#sec:volumes){reference-type="eqref"
reference="sec:volumes"}
and [\[sec:differences\]](#sec:differences){reference-type="eqref"
reference="sec:differences"}---either one---now having

$$\label{eq:multigroup-eigenvalue-discrete}
 R \cdot \boldsymbol{\phi} = \frac{1}{k_\text{eff}} F \cdot \boldsymbol{\phi}$$
where the size of the problem is $N\cdot G$, and the elements of the
eigenvector are $\phi(i,g)$:

$$\boldsymbol{\phi} =
\begin{bmatrix}
 \phi(1,1) \\
 \phi(1,2) \\
 \vdots \\
 \phi(1,G) \\
 \phi(2,1) \\
 \vdots \\
 \phi(2,G) \\
 \phi(3,1) \\
 \vdots \\
 \phi(N,G) \\
\end{bmatrix}$$

Matrices $R$ and $F$ are square $NG \times NG$ matrices. Their elements
are real numbers, namely all the coefficients that multiply each of the
fluxes $\phi(i,g)$ that were derived in
sections [\[sec:volumes\]](#sec:volumes){reference-type="eqref"
reference="sec:volumes"}
and [\[sec:differences\]](#sec:differences){reference-type="eqref"
reference="sec:differences"}. If a coefficient is not divided
by $k_\text{eff}$, it is incorporated into $R$ and into $F$ otherwise.
The row corresponds to the index of the cell and group of the
discretized equation it appears in, and the column corresponds to the
index of the cell and group of the flux that it is multiplying to.

In the discretized
problem [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"}, in addition to some
approximation of the $G$ eigenvalue and eigenvector pairs present in the
continuous
problem [\[eq:multigroup-eigenvalue\]](#eq:multigroup-eigenvalue){reference-type="eqref"
reference="eq:multigroup-eigenvalue"}, there appear $N\cdot(G-1)$ new
eigenvalues. Nonetheless, the three mathematical properties stated above
also hold for
equation [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"}.


The direct storage of a matrix of size $NG \times NG$ for values of $N$
of interest in nuclear reactor analysis is intractable with
computational resources available nowadays[^2]. Luckily, matrices $R$
and $F$ are sparse---while $R^\star$ and $B^\star$ in general are
not---as both the discretized multigruoup diffusion equations and the
boundary conditions depend only on first neighbors in space and on the
fluxes of the other energy groups in the spatial cell evaluated.
However, if the generalized eigenvalue
problem [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"} is to be written as a
standard eigenvalue problem either as

$$F^{-1} R \cdot \boldsymbol{\phi} = \frac{1}{k_\text{eff}} \cdot \boldsymbol{\phi}$$
or as

$$R^{-1} F \cdot \boldsymbol{\phi} = {k_\text{eff}} \cdot \boldsymbol{\phi}$$
in addition to the problem of finding inverse matrices---that, by the
way, may not exists as is the usual case for $F^{-1}$---there appears at
least one big dense matrix that cannot be efficiently handled.
Therefore, to solve the neutron steady-state diffusion equation, methods
for solving a generalized eigenvalue problem in the form of
equation [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"} ought to be used.


Milonga uses the PETSc[^3] [@petsc-user-ref; @petsc-efficient] library
to build and handle $\boldsymbol{\phi}$, $R$ and $F$, and the
SLEPc [@Hernandez:2005:SSF; @Hernandez:2003:SSL] library to solve the
generalized eigenvalue problem given by
equation [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"} to obtain $k_\text{eff}$
and the associated critical reactor flux distribution. SLEPc solves
either the standard or the generalized eigenvalue problem maintaining
the sparsity of the formulation using a variety of efficient iterative
methods, from which milonga's user can choose. In particular, given

$$\label{eq:slepc-generalized}
A \ensuremath\mathbf{x} = \lambda B \ensuremath\mathbf{x}$$ SLEPc tries
to find combinations $(\lambda, \ensuremath\mathbf{x})$ that make the
equation hold. Moreover, it can be told to find only a few large or
small eigenvalues, or to search for solutions located only in a certain
region of the complex plane. As in the diffusion problem, due to the
properties of
equation [\[eq:standard-eigenvalue\]](#eq:standard-eigenvalue){reference-type="eqref"
reference="eq:standard-eigenvalue"} stated in page , only one
eigenvalue-eigenvector pair is the solution sought for in the neutron
diffusion equation, this feature is particularly convenient.


On the one hand, depending on the kind of physical problem being solved,
some solution methods are more suitable than others. On the other hand,
some methods work better---and may even not work at all---when one
of $R$ or $F$ takes place of the matrix $A$ in SLEPc's
equation [\[eq:slepc-generalized\]](#eq:slepc-generalized){reference-type="eqref"
reference="eq:slepc-generalized"}. Thus, essentially two versions of the
generalized problem can be used, either

$$\label{eq:direct_keff}
 F \cdot \boldsymbol{\phi} = k_\text{eff} \, R \cdot \boldsymbol{\phi}$$
or

$$\label{eq:inverse_keff}
 R \cdot \boldsymbol{\phi} = \frac{1}{k_\text{eff}} \, F \cdot \boldsymbol{\phi}$$

In the first case, $\lambda = 1/k_\text{eff}$ and the eigenvalue sought
for is the smaller in magnitude, and conversely.
Equation [\[eq:direct\_keff\]](#eq:direct_keff){reference-type="eqref"
reference="eq:direct_keff"} as the "direct $k_\text{eff}$ formulation"
and
equation [\[eq:inverse\_keff\]](#eq:inverse_keff){reference-type="eqref"
reference="eq:inverse_keff"} is referred to as the "inverse
$k_\text{eff}$ formulation". The user can request to solve the problem
using either one, depending on the particularities of the problem and
the solver method.

### Iterations tolerance

All the eigensolvers provided by SLEPc are iterative in nature, meaning
that the solutions are (usually) improved at each iteration until they
are sufficiently accurate, that is, until convergence is
achieved [@slepc-users-manual]. Thus, a criterion for stopping the
iterative process is needed. It can be given as a maximum number of
iterations or as a tolerance for the solution convergence.

Note that this tolerance is a measure of the error committed by the
numerical method used to solve the eigenvalue
problem [\[eq:inverse\_keff\]](#eq:inverse_keff){reference-type="eqref"
reference="eq:inverse_keff"}
(or [\[eq:direct\_keff\]](#eq:direct_keff){reference-type="eqref"
reference="eq:direct_keff"}). In addition, there exist an error related
to the discretization both in space and energy of the diffusion
equation [\[eq:diffusion-continuous\]](#eq:diffusion-continuous){reference-type="ref"
reference="eq:diffusion-continuous"}, that in turn has uncertainties in
the macroscopic cross-sections. Moreover, the diffusion equation is an
approximation that holds only up to a certain degree only in special
cases of nuclear reactor geometries and material compositions. Keep in
mind
figure [\[fig:uncertainties\]](#fig:uncertainties){reference-type="ref"
reference="fig:uncertainties"} with this respect.


SLEPc considers that the eigenvalue problem is solved whenever a certain
residual is less than some small positive value. There are three
possible criteria, namely absolute residual, relative residual with
respect to the eigenvalue and relative residual with respect to the
matrix norms. As the selection of this tolerance is of vital importance
for the performance of the numerical solution of the diffusion problem,
a trade-off between high accuracy and a reasonable computation speed is
needed.

According to the SLEPc manual [@slepc-users-manual], the provided
convergence criteria for the generalized eigenvalue
problem [\[eq:slepc-generalized\]](#eq:slepc-generalized){reference-type="eqref"
reference="eq:slepc-generalized"} are based on the absolute residual
vector defined as

$$r = A \cdot \tilde{\ensuremath\mathbf{x}} - \tilde{\lambda} \, B \cdot \tilde{\ensuremath\mathbf{x}}$$
where $\tilde{\lambda}$ and $\tilde{\ensuremath\mathbf{x}}$ represent
the computed eigenpair. For Hermitian problems---which is not the case
for $R$ and $F$---it is possible to prove, as stated
in [@slepc-users-manual], that

$$\left| \lambda - \tilde{\lambda} \right| \leq \| r \|_2$$

This is not the case in general problems, so there is no simple
relationship between the residual norm and the actual error committed by
the numerical method. Nevertheless, $\| r \|_2$ stills should give a
measure of the distance between the real eigenvalue $\lambda$ and the
approximate solution $\lambda$.


As $\| r \|$ depends on the size of the problem, in general a finer
discretization will lead to more iterations for the eigenvalue problem.
To avoid this escalation, milonga uses as a measure of the convergence
the residual norm with respect to the matrices norm, namely

$$\label{eq:eigen-converged}
 \| r_\text{rel} \| = \frac{\| r \|_2}{\| A \|_2 + |\lambda| \| B \|_2 }$$
so convergence is considered as achieved whenever $\| r_\text{rel} \|$
is less than a certain tolerance. The default tolerance value used by
milonga is the default used by SLEPc that is $10^{-7}$. The impact of
this figure on the overall behavior with respect of accuracy and running
time of milonga when solving a particular problem has to be determined
by the user in a by-case basis.

### Flux and power {#sec:power}

In the absence of external independent sources, the linear diffusion
problem [\[eq:diffusion-continuous\]](#eq:diffusion-continuous){reference-type="eqref"
reference="eq:diffusion-continuous"} is homogeneous, and so is the
associated discretized eigenvalue
problem [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"}. Therefore the solution
is defined up to a multiplicative constant. That is, if
$\boldsymbol{\phi}_1$ is a solution, then
$\boldsymbol{\phi}_2 = \alpha \boldsymbol{\phi}_1$ is also a solution.
To be able to compare solutions obtained by different methods, a
normalization procedure is desired. Milonga computes the
solution $\boldsymbol{\phi}$ of the eigenvalue
problem [\[eq:multigroup-eigenvalue-discrete\]](#eq:multigroup-eigenvalue-discrete){reference-type="eqref"
reference="eq:multigroup-eigenvalue-discrete"} such that the mean value
is equal to one unit of flux, i.e.

$$\frac{\displaystyle \int_V \int_0^{\infty} \phi(\ensuremath\mathbf{r},E) \, dE \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_V d^m \ensuremath\mathbf{r}} = 1$$

Using
equations [\[eq:energy-integrated-flux\]](#eq:energy-integrated-flux){reference-type="eqref"
reference="eq:energy-integrated-flux"}
and [\[eq:flux-in-volumes\]](#eq:flux-in-volumes){reference-type="eqref"
reference="eq:flux-in-volumes"}, this condition can be written as

$$\begin{aligned}
 \frac{\displaystyle \int_V \sum_{g=1}^G \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_V d^m \ensuremath\mathbf{r}} &= 1 \\
 \frac{\displaystyle \sum_{g=1}^G \int_V  \phi(\ensuremath\mathbf{r}, g) \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_V d^m \ensuremath\mathbf{r}} &= 1 \\
 \frac{\displaystyle \sum_{g=1}^G \sum_{i=1}^N \phi(i, g) \int_{V_i} \, d^m\ensuremath\mathbf{r}}{\displaystyle \int_V d^m \ensuremath\mathbf{r}} &= 1\end{aligned}$$
that for uniformly distributed cells, reduces to

$$\begin{aligned}
 \sum_{g=1}^G \sum_{i=1}^N \phi(i, g) &= 1\end{aligned}$$

In real reactors, feedback effects---such as temperature or
xenon---render the problem nonlinear and the solution cannot be
modulated at will, but the flux distribution is unique. To be able to
cope with these effects, a desired power level has to be given in order
to normalize the flux distribution to the correct value. If the power
setpoint is $P^\star$, for a three-dimensional reactor, this condition
reads

$$\int_V \int_0^{\infty} E\Sigma_f(\ensuremath\mathbf{r}, E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE \, d^3\ensuremath\mathbf{r} = P^\star$$
where $E\Sigma_f$ is the product of energy released in one single
fission times the macroscopic fission cross section. Note that if the
desired power level $P^\star$ refers to the prompt fission power, then
$E\Sigma_f$ has to include the prompt energy only, while if $P^\star$ is
the total power---fission plus delayed---then all the contributions to
the fission energy by the fission product decay chain should be taken
into account. Otherwise, a correction factor to take into account
delayed power has to be incorporated either in $P^\star$ or
in $E\Sigma_f$.


Three-dimensional reactors occupy a finite volume in space, and thus
$P^\star$ is the reactor power, either total or fission according to the
definition of $E\Sigma_f$. But one and two-dimensional reactors are
essentially infinite, and thus the dissipated power is also infinite.
For one dimension, $P^\star$ is the power density per unit area
perpendicular to the $x$ axis, and for two dimensions it is the power
density per unit length perpendicular to the $x$-$y$ plane. Note that
$E\Sigma_f \cdot \phi$ has units of energy $\cdot$ length$^{-3}$. Units
consistency is expected.


When a power setpoint is defined, the flux is normalized such that

$$\begin{aligned}
\int_V \int_0^{\infty} E\Sigma_f(\ensuremath\mathbf{r}, E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE \, d^3\ensuremath\mathbf{r} &= P^\star \\
\int_V \sum_{g=1}^{G} E\Sigma_f(\ensuremath\mathbf{r}, g) \cdot \phi(\ensuremath\mathbf{r}, g)  \, d^3\ensuremath\mathbf{r} &= P^\star \\
\sum_{g=1}^{G} \sum_{i=1}^N \int_{V_i} E\Sigma_f(\ensuremath\mathbf{r}, g) \cdot \phi(\ensuremath\mathbf{r}, g)  \, d^3\ensuremath\mathbf{r} = P^\star \\
\sum_{g=1}^{G} \sum_{i=1}^N E\Sigma_f(i, g) \cdot \phi(i, g)  \int_{V_i} \, d^3\ensuremath\mathbf{r} = P^\star \\\end{aligned}$$
where
equation [\[eq:flux-in-volumes\]](#eq:flux-in-volumes){reference-type="eqref"
reference="eq:flux-in-volumes"} was used, and the multigroup
property $E\Sigma_f$ was defined---as with $\Sigma_t$ in
equation [\[eq:sigma\_t\_weighted\]](#eq:sigma_t_weighted){reference-type="eqref"
reference="eq:sigma_t_weighted"}---as

$$E\Sigma_f(\ensuremath\mathbf{r},g) = \frac{\displaystyle  \int_{E_g}^{E_{g-1}} E\Sigma_f(\ensuremath\mathbf{r},E) \cdot \phi(\ensuremath\mathbf{r},E) \, dE}{\displaystyle \int_{E_g}^{E_{g-1}} \phi(\ensuremath\mathbf{r},E) \, dE}$$
and then associated to a particular cell $E\Sigma_f(i,g)$ with one of
the methods discussed in
section [2.2](#sec:xs-association){reference-type="ref"
reference="sec:xs-association"}.


Note that the integrals over the spatial domain were replaced by sums by
using
equation [\[eq:flux-in-volumes\]](#eq:flux-in-volumes){reference-type="eqref"
reference="eq:flux-in-volumes"}, and thus are mathematically exact from
this point of view. However, if using finite differences, $\phi(i,g)$ is
not the mean value of the flux in a cell but the value of the flux at
the cell center. Therefore, in this case the replacement of the integral
by the sum is exact up to order $\Delta x^2$---i.e. as accurate as an
integration by the trapezoids method---given by the the
approximation [\[eq:approx-meanflux-centerflux\]](#eq:approx-meanflux-centerflux){reference-type="eqref"
reference="eq:approx-meanflux-centerflux"}.

Xenon poisoning
---------------

Depending on the flux level distribution, fission power reactors get
poisoned by the appearance of fission products. Of all the possible
isotopes, the mos important from the neutronic point of view
is $^{135}$Xe, and it is usually the only one taken into account in core
calculations. Milonga is capable of handling the effects of xenon in the
power distribution, using the models and equations discussed below.

The mechanisms of appearance and removal of the isotope $^{135}$Xe are
illustrated in figure [\[fig:xenon\]](#fig:xenon){reference-type="ref"
reference="fig:xenon"}. On the one hand, it can appear as a direct
fission product with a yield $\gamma_\text{Xe}$---the mean number of
nuclides of this isotope that appear promptly after a single
fission---or as a daughter of $^{135}$I, that decays into xenon with a
time constant $\lambda_\text{I}$. In turn, $^{135}$I appears as a
daughter of $^{135}$Te that is a fission product. However, as the beta
decay of tellurium is so fast with respect to the decay of iodine,
$^{135}$I is considered to appear as a prompt fission product. On the
other hand, $^{135}$Xe disappears because of beta decay into caesium and
because of neutron absorption, transforming into $^{136}$Xe at a rate
that depends on the local flux.

![[\[fig:xenon\]]{#fig:xenon label="fig:xenon"}Mechanisms of appearance
and removal of $^{135}$ as modeled by milonga.](xenon)

At a point $\ensuremath\mathbf{r}$ in the reactor, the iodine
concentration $\text{I}(\ensuremath\mathbf{r}, t)$ is assumed to follow
the balance equation

$$\label{eq:iodine-continuous}
 \frac{\partial \text{I}(\ensuremath\mathbf{r}, t)}{\partial t} = \gamma_\text{I} \cdot \int_0^{\infty} \Sigma_f(\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE  - \lambda_\text{I} \cdot \text{I}(\ensuremath\mathbf{r},t)$$
and the $^{135}$Xe is given by

$$\begin{aligned}
\label{eq:xenon-continuous}
 \frac{\partial \text{Xe}(\ensuremath\mathbf{r}, t)}{\partial t} =& \,\gamma_\text{Xe} \cdot \int_0^{\infty} \Sigma_f(\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE + \lambda_\text{I} \cdot \text{I}(\ensuremath\mathbf{r},t) \nonumber \\
& \quad\quad - \lambda_\text{Xe} \cdot \text{Xe}(\ensuremath\mathbf{r},t) - \int_0^{\infty} \Sigma_{a\text{Xe}} (\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE\end{aligned}$$

The total fissions can be computed from the group fluxes as

$$\int_0^{\infty} \Sigma_f(\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE = \frac{1}{\nu} \sum_{g=1}^{G} \nu\Sigma_f(\ensuremath\mathbf{r}, g, t) \cdot \phi(\ensuremath\mathbf{r}, g, t)$$
where a mean number of neutrons emitted per fission was used in order to
avoid introducing a new nuclear parameter $\Sigma_f$.

The macroscopic cross section of xenon can be separated into

$$\Sigma_{a\text{Xe}}(\ensuremath\mathbf{r}, E, t) = \sigma_{a\text{Xe}}(E) \cdot \text{Xe}(\ensuremath\mathbf{r}, t)$$
as the xenon is assumed to be so diluted that its microscopic absorption
cross section $\sigma_{a\text{Xe}}$ depends only on energy. Therefore,
the rate of removal of $^{135}$Xe due to neutron absorption in
equation [\[eq:xenon-continuous\]](#eq:xenon-continuous){reference-type="eqref"
reference="eq:xenon-continuous"} can be written as

$$\begin{aligned}
 \int_0^{\infty} \Sigma_{a\text{Xe}} (\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE &= \int_0^{\infty} \sigma_{a\text{Xe}}(E) \cdot \text{Xe}(\ensuremath\mathbf{r}, t) \cdot (\ensuremath\mathbf{r}, E, t) \cdot \phi(\ensuremath\mathbf{r}, E, t) \, dE \\
&= \text{Xe}(\ensuremath\mathbf{r}, t) \cdot \sum_{i=1}^{G} \sigma_{a\text{Xe}}(g) \cdot \phi(\ensuremath\mathbf{r},g,t)\end{aligned}$$

The fixed point of
equations [\[eq:iodine-continuous\]](#eq:iodine-continuous){reference-type="eqref"
reference="eq:iodine-continuous"}
and [\[eq:xenon-continuous\]](#eq:xenon-continuous){reference-type="eqref"
reference="eq:xenon-continuous"} gives the steady state xenon
concentration $\text{Xe}^\star(\ensuremath\mathbf{r})$, that in terms of
the multigroup fluxes is

$$\text{Xe}^\star(\ensuremath\mathbf{r}) = \frac{\displaystyle \big(\gamma_\text{I} + \gamma_\text{Xe}\big) \cdot  \frac{1}{\nu} \sum_{g=1}^G \nu\Sigma_f(\ensuremath\mathbf{r},g) \cdot \phi(\ensuremath\mathbf{r},g) }{\displaystyle \lambda_\text{Xe} + \sum_{g=1}^G \sigma_{a\text{Xe}}(g) \cdot \phi(\ensuremath\mathbf{r},g)}$$

In particular, at the $i$-th cell center the steady-state xenon
concentration is

$$\label{eq:xenon}
 \text{Xe}^\star(i) = \frac{\displaystyle \big(\gamma_\text{I} + \gamma_\text{Xe}\big) \cdot  \frac{1}{\nu} \sum_{g=1}^G \nu\Sigma_f(i, g) \cdot \phi(i,g) }{\displaystyle \lambda_\text{Xe} + \sum_{g=1}^G \sigma_{a\text{Xe}}(g) \cdot \phi(i,g)}$$

This concentration distribution usually modifies the macroscopic cross
sections that appear as coefficients in the diffusion
equation [\[eq:diffusion-continuous-P\]](#eq:diffusion-continuous-P){reference-type="eqref"
reference="eq:diffusion-continuous-P"}, rendering the problem
non-linear. Thus, a iterative procedure is needed to solve for the
actual flux distribution taking into account the effects of xenon
poisoning, as discussed in page . Note that the xenon concentration
depends on the absolute value of the neutron flux, thus a power setpoint
is needed in order to correctly normalize the solution of the eigenvalue
problem to dimensional quantities, as explained in
section [2.4.2](#sec:power){reference-type="ref" reference="sec:power"}.
Concentrations have units of inverse cubed length, decay constants
inverse time and microscopic cross sections inverse squared length.
Yields and the mean number of neutrons emitted by a single fission are
dimensionless numbers.

\bibliographystyle{unsrt}
Input preparation {#cap:input}
=================


\hspace{\fill}
\sf 
\footnotesize
Science is what we understand well enough to explain to a computer.\
Art is everything else we do.

*Donald Knuth, Foreword to the book "A=B", 1996*


This chapter is about imagination. It is about asking a dumb computer to
solve a set of equations that somehow represent a nuclear reactor. And
according to milonga's design basis, chances are that this nuclear
reactor exists nowhere but in your mind. If after reading this chapter
you still cannot ask milonga to help you with your fancy
state-of-the-art reactor design, please do one of the following two
things: contact the author to comply about the lack of the features you
feel you are needing, or---even better---hack into the source code and
add the feature by yourself. And of course, distribute the modified code
under the terms of the GNU Public License or---even better---share your
modifications with the author to have them appended into the official
milonga distribution.


Although milonga has a lot of options that can be tweaked, it comes with
a great deal of default values that should satisfy most needs. That is
why input files can be either very small or arbitrary large. For
instance, the following input is the very smallest one that gets the
program to perform a useful calculation:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PROBLEM DIMENSIONS 1 GROUPS 1
    x_bare_length = 100
    x_cells = 100
    ZONE fuel MATERIAL fuel
    MATERIAL fuel D_1 1 SigmaA_1 1e-2 nuSigmaF_1 1.2e-2

The rest of this chapter is devoted to explaining how to build your own
input files to, hopefully, solve the problem you have at hand.

Input files {#sec:input}
-----------

Being built on top of wasora's framework, milonga provides all the basic
functionality in wasora plus specific options and settings for milonga.
A detailed description of the keywords associated to the common
framework can be found in wasora's documentation. Nevertheless, a quick
explanation of the basic wasora's features is given here for
completeness.

Milonga's normal behavior is to read one or more input files defining
the problem, solve it and then output the requested results---if
any---in different possible ways. The input files are plain text files
containing keywords, parameters and or numerical data. There should be
one main input file than can optionally include further input files.

There are primary and secondary keywords. Primary keywords should appear
as the first token in a line. Secondary keywords are usually optional
and should appear in the same line after a primary keyword. One primary
keyword may be followed by zero or more secondary keywords. The order of
the secondary keywords does not matter, while the order in which the
primary keywords appear in the input file may or may not matter.
Keywords are case insensitive, but for clarity they are written in
uppercase throughout this document, except for the secondary keywords of
`MATERIAL` that define the nuclear properties.

Both primary and secondary keywords may take zero or more parameters.
Each parameter should be a single token delimited by spaces or tabs. If
spaces are needed as a part of the parameter, it should be enclosed in
double quotes.

Blank lines and characters appearing after a hash '\#' symbol are
treated as separators or comments and are ignored. For readability
reasons, if a definition is too long to fit neatly in a single line, it
can be enclosed in braces '{' and '}' spanning different file lines.
From milonga's point of view, whatever appears enclosed in braces is
treated as belonging to the same line.


The following examples illustrates the general keyword syntax:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # this is an example of the general syntax expected in input files
    PRIMARY_KEYWORD KEYWORD1 argument1 KEYWORD2 "one single argument"

    # the above line could also be written as
    PRIMARY_KEYWORD {
      KEYWORD1 argument1       # with comments between lines
      KEYWORD2 "one single argument"
    }


There are keywords that are common to all the packages belonging to the
wasora suite---i.e. the common framework---and keywords that are
specific to milonga. Special stress is given to specific keywords,
although as some reference to common keywords is done, a basic
explanation about the framework's features is also provided.


Some information is given to milonga by giving arguments to keywords,
and others by defining special variables to take certain values. This is
the case for information about the problem that may either be given from
external sources or needed after the computation is given. For example,
the number of energy groups is a rather fixed parameter and is entered
as a keyword argument, while the number of spacial cells in a certain
direction is given as a variable because this value may be part of a
parametric study and/or needed to be used as an output or to request the
flux distribution at the computed points to avoid interpolation.


In any case, most of the parameters are entered by using algebraic
expressions accepted by the wasora framework. These expressions are a
combination of numeric constants, variables, functions and operators. A
complete description can be found in wasora's documentation, but valid
expressions are fairly intuitive:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    x_bare_length = 100
    x_cells = round(random(1, 100))
    pi = 4 * atan(1)
    buckling = (pi/x_bare_length)^2

    PRINT TEXT "buckling = " buckling
    PRINT TEXT "delta x =  " x_bare_length/x_cells

In this example, `x_bare_length` and `x_cells` are special variables
that are used by milonga to define the problem, as explained below. The
bare length is fixed to one hundred, but the number of cells is a random
integer number between one and one hundred. Then, the number $\pi$ is
computed by using the arc-tangent of one, and finally the geometric
buckling is calculated. As an example of information given as keywords
arguments, consider the last two lines. Keyword `PRINT` accepts
expressions as arguments. The first instruction has a single variable as
the argument to be printed, but the second one has a quotient.

Note that the expressions that appear after the equal sign '=' in
assignments can have blanks and spaces as desired, as long as the whole
expression fits in a single line (or is enclosed in brackets). When
appearing as keyword arguments, there should be no intermediate blanks
and the whole expression should be a single token. If not, the
expression should be enclosed in double quotes.


In any location of an input file, another input file can be included
with the keyword


    INCLUDE file


where `file` is either a relative or absolute path. The effect is as if
the included file was pasted right in the place where the keyword was
inserted in the original input file. Nested inclusions, i.e. `INCLUDE`
keywords in included files, are possible.


During this chapter, optional constructions are presented as enclosed in
square brackets '\[' and '\]'. Whenever only one of several keywords is
to be selected, all the possible keywords are enclosed in
parenthesis '(' and ')' separated by a vertical line '$|$'. Keywords
arguments are written in lowercase, and if the content and meaning is
not explained in the surrounding paragraphs, it should be inferred by
the context.


As with the example about variables and keyword arguments above, the
rest of this chapter contains partial examples of the usage of the
feature being discussed. However, they are just illustrative and by no
means provide real or ready-to-use inputs.
Chapter [4](#cap:examples){reference-type="ref"
reference="cap:examples"} contains several stand-alone examples of input
files that span a lot of different kinds of problems and conditions.
They are arranged in increasing order of complexity, and the inputs
shown should be fairly intuitive and self explanatory. If there is a
doubt about the usage of a certain feature, these complete examples
should give the correct answer about how to use the code.

Problem definition
------------------

The first definition that ought to be provided to milonga is the number
of spatial dimensions and the number of energy groups. This information
is given by the primary keyword `PROBLEM`, that takes the secondary
keywords `DIMENSIONS` and `GROUPS`. Each one takes one integer argument,
thus defining the dimension of the problem phase space


    PROBLEM DIMENSIONS m GROUPS g


This information is mandatory and should appear before any other of
milonga's specific keywords. For example, to define a two-dimensional
problem with three energy groups:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PROBLEM DIMENSIONS 2 GROUPS 3


To define the geometry, you start with an $m$-dimensional bare prism and
then you define one or more material zones inside it. This prism is
characterized by $m$ lengths, called the bare lengths. These lengths are
entered as three variables, one for each direction:


    x_bare_length
    y_bare_length
    z_bare_length


For two-dimensional problems, `z_bare_length` is ignored and does not
need to be entered. The same applies for `y_bare_length` in
one-dimensional problems. The units of these variables are to be decided
by the user. However, it is important to note that the entered values
ought to have the same units, consistent with the cross sections entered
in the material definition
(section [3.5](#sec:materials){reference-type="ref"
reference="sec:materials"}). That is, if the macroscopic cross sections
are in cm$^{-1}$, then the bare lengths should be in cm. For
non-dimensional problems, the consistency between the
nondimensionalization of the lengths and the cross sections is up to the
user.


To define a three-dimensional bare prism from which to start defining
the reactor geometry with a size of 100$\times$100$\times$100 of a
certain units consistent with the cross sections, say centimeters, enter

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize,frame=single}
    x_bare_length = 100
    y_bare_length = 100
    z_bare_length = 100

The expressions can be any valid mathematical expression. For example,
in a two-dimensional problem, to make sure the bare prism is a always a
square, write

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize,frame=single}
    x_bare_length = a
    y_bare_length = x_bare_length

The last thing that milonga should know about the problem is the spatial
nodalization. This version only supports uniformly distributed cells, so
the number of cells in each direction has to be given. To do so, use the
following variables:


    x_cells
    y_cells
    z_cells


These variables should be non-dimensional integers. If any of these
values is non-integral, it will be rounded to the nearest integer. A
valid construction is

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    x_cells = 100
    y_cells = 50
    z_cells = 40

Again, any expression can be used, so the actual nodalization can be
entered as a cell length instead of cell numbers as:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    delta_x = 1
    delta_y = 2
    delta_z = 2.5
    x_cells = x_bare_length/delta_x
    y_cells = x_bare_length/delta_y
    z_cells = y_bare_length/delta_z

This information should be entered before the definition of the material
zones, that is explained in
section [3.3](#sec:zones){reference-type="ref" reference="sec:zones"}.


Milonga assumes that there is a single fissionable isotope or,
equivalently, that all the available isotopes have the same fission
spectrum. The default fission spectrum is thermal, i.e., all the
neutrons are born in the last energy group. If this behavior is to be
changed, the elements of the vector

    chi

have to be entered, according to
equation [\[eq:chi\_weighted\]](#eq:chi_weighted){reference-type="eqref"
reference="eq:chi_weighted"}. Elements of a vector in wasora are
accessed by appending an underscore '\_' plus the element number,
starting in one. Thus, in a three group problem, the default fission
spectrum is taken as if the user had entered

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    chi_1 = 0
    chi_2 = 0
    chi_3 = 1

It is the user's responsibility to ensure that the sum of the components
of vector $\ensuremath\mathbf{\chi}$ to sum up to one, as required by
the normalization condition given in
equation [\[eq:chi\_normalized\]](#eq:chi_normalized){reference-type="eqref"
reference="eq:chi_normalized"}.

Zones {#sec:zones}
-----

The actual reactor geometry is defined by entering one or more material
zones. This is accomplished by using the `ZONE` primary keyword. It
takes one mandatory argument that defines an identification string and
the mandatory `MATERIAL` secondary keyword that associates a material to
the zone. Then, there may be zero or more optional secondary keywords
with arguments that define the geometric place of the zone. There is
also one optional secondary keyword `INCREMENTAL`.


    ZONE name MATERIAL mat [X_MIN expr] [X_MAX expr] [Y_MIN expr] [Y_MIN expr] [Z_MIN expr] [Z_MAX expr] [X_CENTER expr] [Y_CENTER expr] [INNER_RADIUS expr] [OUTER_RADIUS expr] [INCREMENTAL]


The identification `name` should be unique for each zone. Milonga gives
an error if two zones have the same name. The material should be one
defined with a `MATERIAL` primary keyword
(section [3.5](#sec:materials){reference-type="ref"
reference="sec:materials"}), and the name is case sensitive. It is not
needed to have the material defined before the zone in the input file,
so the material definitions can be left to the end of the file---or
maintained in an `INCLUDE`'d file---for readability.


By default, if no geometry keywords are entered, the zone spans the
whole bare prism. Currently, only prisms with faces parallel to the bare
prisms and circles in the $x$-$y$ plane are supported in this version.
The limits can be given as expressions, as long as there are no
intermediate blanks. Use quotes if needed, as explained in
section [3.1](#sec:input){reference-type="ref" reference="sec:input"}.
The units should be the same as the units of `x_bare_length`.

If the min or max keywords are given, the zone is limited by a plane
perpendicular to the specified axis, otherwise the zone extends up to
the bare prism limit in that direction. In two dimensions, a circle can
be entered by giving the coordinates of the center and the radius. An
annulus is entered by giving also a value for the inner radius.

Zones entered last in the input file overwrite zones previously defined.
If the optional keyword `INCREMENTAL` is given, then the cross sections
of the material associated are algebraically added to whatever zones are
existing under the new one.


For example, to define one single region that spans the whole bare prism
with a material named "fuel", write

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    ZONE core MATERIAL fuel

To define a 20$\times$10 (in whatever units `x_bare_length` is)
rectangle centered in the bare square, enter

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    ZONE rec MATERIAL fuel { X_MIN x_bare_length/2-10 X_MAX x_bare_length/2+10
                             Y_MIN y_bare_length/2-5  Y_MAX y_bare_length/2+5 }

If a zone is to have a circular shape, `X_CENTER`, `Y_CENTER` and
`OUTER_RADIUS` have to be given. To make an annulus, provide also
`INNER_RADIUS`. To illustrate this and the overlapping logic, consider
these two equivalent definitions:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    ZONE refl MATERIAL reflector X_CENTER 200 Y_CENTER 200 OUTER_RADIUS 200
    ZONE core MATERIAL fuel      X_CENTER 200 Y_CENTER 200 OUTER_RADIUS 180

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    ZONE core MATERIAL fuel      X_CENTER 200 Y_CENTER 200 OUTER_RADIUS 180
    ZONE refl MATERIAL reflector X_CENTER 200 Y_CENTER 200 OUTER_RADIUS 200 INNER_RADIUS 180

Note that currently milonga does not support internal vacuum zones, i.e.
zones that are surrounded with materials but that are not assigned to
any material themselves.

Boundary conditions {#boundary-conditions}
-------------------

To set the problem boundary conditions, milonga detects the external
surfaces of the zones defined---that may or may not coincide with the
bare prism external surfaces---and applies the equations described in
section section [2.3](#sec:boundary){reference-type="ref"
reference="sec:boundary"} to the corresponding cells. The default
boundary condition is to require null flux at the external surfaces.
This version supports either this default behavior or mirror conditions
on planar surfaces. Circular external surfaces always have null flux
conditions.

To change from null to mirror conditions, use the primary keyword
`BOUNDARY_CONDITIONS`. This keyword takes optional secondary keywords
that define the kind of boundary condition for each axis and for each
direction


    BOUNDARY_CONDITIONS [X_MIN (NULL | MIRROR)] [X_MAX (NULL | MIRROR)] [Y_MIN (NULL | MIRROR)] [Y_MAX (NULL | MIRROR)] [Z_MIN (NULL | MIRROR)] [Z_MAX (NULL | MIRROR)]


The order in the input file at which the `BOUNDARY_CONDITIONS` keyword
appears is irrelevant. If one secondary keyword is omitted, the
corresponding boundary condition is taken as `NULL` by default. For
example, in one dimension to have mirror conditions in the left and null
flux in the right, use

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    BOUNDARY_CONDITIONS X_MIN MIRROR

Even though external circular surfaces are forced to have null flux, the
bare prism surfaces may belong to a circle and have mirror conditions as
well. For example, to implement a quarter of a symmetric circular core
define

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    PROBLEM DIMENSIONS 2 GROUPS 1
    x_bare_length = 50
    y_bare_length = 50

    ZONE core MATERIAL fuel X_CENTER 0 Y_CENTER 0 OUTER_RADIUS 50

    BOUNDARY_CONDITIONS X_MIN MIRROR Y_MIN MIRROR

Keep in mind that the requirement of null fluxes in a discrete domain
may lead to artificial negative cell fluxes depending on the relative
position of the cell center and the external surface, as explained in
page  and depicted in
figure [\[fig:bc-planarneg\]](#fig:bc-planarneg){reference-type="ref"
reference="fig:bc-planarneg"}.

Materials {#sec:materials}
---------

The zones discussed in section [3.3](#sec:zones){reference-type="ref"
reference="sec:zones"} have one associated material, that is defined
with the `MATERIAL` primary keyword. Note that one zone has only one
associated material, but one material can be associated to more than one
zone. The order at which the materials are defined in the input file is
irrelevant.

    MATERIAL name [D_1 expr] [D_2 expr] [...] ([SigmaA_1 expr] [SigmaA_2 expr] [...] | [SigmaT_1 expr] [SigmaT_2 expr] [...]) [nuSigmaF_1 expr] [nuSigmaF_1 expr] [...] [SigmaS_1.1 expr] [SigmaS_1.2 expr] [...] [SigmaS_2.1 expr] [SigmaS_2.2 expr] [...] [...] [ESigmaF_1 expr] [ESigmaF_2 expr] [...]

The identification `name` should be unique for each material. Milonga
gives an error if two materials have the same name. The name is case
sensitive. Secondary keywords define the nuclear parameters that
characterize the material. If no information is entered about a certain
parameter, it is taken as being identically equal to zero.

These secondary keywords have a main name that represents the kind of
nuclear parameter, followed by an underscore '\_' and then a group
number, except `SigmaS` that have two group numbers separated by a point
'.'. Table [\[tab:xs\]](#tab:xs){reference-type="ref"
reference="tab:xs"} lists the main names, the description and the units
of the expected argument. The dots '\[...\]' in the above syntax mean
that one keyword for each energy group is expected. For the scattering
cross sections, all the $G^2$ groups combinations can be entered.

    Keyword        Symbol      Description                                                                                                        Units
  ------------ --------------- ------------------------------------------------------------------------------------------------------ ------------------------------
      `D`            $D$       Diffusion coefficient                                                                                              length
    `SigmaA`     $\Sigma_a$    Absorption macroscopic cross section                                                                           length$^{-1}$
    `SigmaT`     $\Sigma_t$    Total macroscopic cross section                                                                                length$^{-1}$
   `nuSigmaF`   $\nu\Sigma_f$  Product of the mean number of neutrons born in a fission times the fission macroscopic cross section           length$^{-1}$
    `SigmaS`     $\Sigma_s$    Scattering macroscopic cross section                                                                           length$^{-1}$
   `ESigmaF`     $E\Sigma_f$   Product of the mean energy released in a fission times the fission macroscopic cross section            energy $\cdot$ length$^{-1}$

  : [\[tab:xs\]]{#tab:xs label="tab:xs"}Secondary keywords for the
  primary keyword `MATERIAL`. The symbols are the ones used throughout
  chapter [2](#cap:equations){reference-type="ref"
  reference="cap:equations"}. Units of length have to be consistent with
  the geometry definition and units of energy have to be consistent with
  the power setpoint, if given. Also, the character of the energy
  (prompt or total) should also be consistent with the character of the
  power setpoint.


Milonga uses the total macroscopic cross sections for the multigroup
formulation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"}. However, lattice codes usually
output the absorption cross section instead. As these two parameters are
related by

$$\Sigma_t(g) = \Sigma_a(g) + \sum_{g^{\prime}=1}^G \Sigma_s(g \rightarrow g^\prime)$$
the total cross section can be computed by milonga from the absorption
and scattering cross sections.


The nuclear parameters entered should be consistent with the
mathematical development discussed in
section [2.1.1](#sec:multigroup){reference-type="ref"
reference="sec:multigroup"}. Again, it is stressed that it is the user's
responsibility to assure that the macroscopic cross sections are as
expected by milonga according to its multigroup formulation.


The scattering cross sections keywords are formed with two group numbers
separated by a point. The first number is the energy group of the
impinging neutrons and the second is the energy of the scattered
neutrons. That is to say, `SigmaS_1.2` is the scattering cross section
from group one to group two, i.e., the down-scattering cross section
when taking as usual the first group as the thermal one in a two-group
formulation. Note that when entering the absorption cross section
instead of the total cross section, the value of the self-scattering
cross section is irrelevant, as it appears both in a positive and in a
negative term in
equation [\[eq:diffusion\_multigroup\]](#eq:diffusion_multigroup){reference-type="eqref"
reference="eq:diffusion_multigroup"} and cancels out.


When giving a power setpoint, the `ESigmaF` keywords should also be
given for at least one material. This values represent the product of
the energy released in a single fission times the macroscopic fission
cross section for each group, and is used to compute the reactor power
as explained in section [2.4.2](#sec:power){reference-type="ref"
reference="sec:power"}. It is important to note that consistency is
expected both in the units and in the sense of prompt or total character
of the fission energy and the power setpoint values as discussed in
page .


Of course, nuclear properties may change within the same material
because the parameters on which the cross sections depend may have
non-trivial spatial distributions. This dependence can be taken into
account by employing mathematical expressions instead of numerical
constants as arguments to the secondary keywords of `MATERIAL`. These
expressions can be explicit or implicit functions of the special
variables `x`, `y` and `z`, that will be replaced by milonga with the
components of the position vector $\ensuremath\mathbf{r}$ when
evaluating the cell nuclear properties as explained in
section [2.2](#sec:xs-association){reference-type="ref"
reference="sec:xs-association"}. In general, macroscopic cross sections
depend on certain parameters---such as temperatures or burn-up---that,
in turn, have a certain spatial distribution.

For example, if $\Sigma_a$ depends on the fuel
temperature $T_\text{cool}$ and the coolant temperature $T_\text{cool}$
as a two-dimensional scalar field $f(T_\text{fuel}, T_\text{cool})$ and,
in turn, the spatial distribution of the fuel and coolant temperature is
given by $g(x,y,z)$ and $h(x,y,z)$ repectively, then

$$\Sigma_a(x,y,z) = f\Big(g\left(x,y,z\right),h\left(x,y,z\right)\Big)$$

To implement this in milonga, assume you have defined multidimensional
functions $f(T_\text{fuel}, T_\text{cool})$, $g(x,y,z)$ and $h(x,y,z)$
that interpolate scattered data provided in a text file using wasora's
common framework keywords:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    FUNCTION f(tfuel,tcool) FILE sigmaa.dat
    FUNCTION g(x,y,z)       FILE tfuel-dist.dat
    FUNCTION h(x,y,z)       FILE tcool-dist.dat

Then, the appropriate expression for the absorption cross section within
the material definition would be

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    MATERIAL fuel SigmaA_1 f(g(x,y,z),h(x,y,z)) [...]

One big deal of milonga's flexibility discussed in
section [1.2](#sec:basis){reference-type="ref" reference="sec:basis"}
appears around the collection of available ways to enter this
distributions. The expressions for the nuclear properties should provide
functions of $\ensuremath\mathbf{r}$---through $x$, $y$ and $z$---that
are continuous in the sense that they should be prone to be evaluated at
any arbitrary point of the spatial domain, and not be limited to be
computed only at discrete points. These expressions can use any of the
wasora's features related to function manipulation. Functions can be
provided algebraically or read from scattered data from files or from
external codes trough shared memory objects. Many interpolation methods
for one and more dimensions are available, and further operations on
functions such as integration or differentiation can be done.


To illustrate the joint usage of the `MATERIAL` keyword with continuous
functions, consider the following definition, where a one-dimensional
fuel and temperature distributions are assumed and the cross sections of
the fuel are entered as central values plus partial derivatives
multiplied by parameters increments:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    FUNCTION Tfuel(x) = 400 + 200*sin(x*pi/300)
    FUNCTION Tcool(x) = 250 - 25*cos(x*pi/300)

    MATERIAL fuel {
      D_1        1
      SigmaA_1   "1.0e-2 + 1e-4*(Tfuel(x) - 400) - 5e-5*(Tcool(x) - 250)"
      nuSigmaF_1 "1.2e-2 + 5e-5*(Tfuel(x) - 400)"
    }

Further and complete examples---including multidimensional interpolation
cases---can be found in chapter [4](#cap:examples){reference-type="ref"
reference="cap:examples"}.

Collecting results
------------------

Milonga uses special variables both for gathering information to solve
the problem and for giving back the results to the user. The most
important variable is called


    keff


that contains the effective multiplication factor as obtained by solving
the eigenvalue problem. This variable is defined when starting the
execution and initialized to one, and its content is "filled" with the
correct value after the first `ZONE` keyword. That is to say, any
algebraic expression involving `keff` should appear after the last zone
was defined. For example, in

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    rho1 = (keff-1)/keff
    ZONE fuel MATERIAL fuel
    rho2 = (keff-1)/keff

only variable `rho2` contains the correct reactivity of the problem
being solved, while `rho1` is equal to zero.


The other output variables given by milonga are


    keff_error
    build_time
    solve_time


The first one, `keff_error`, gives an estimation of the absolute error
committed in the solution of the eigenvalue problem, `build_time` that
gives the time needed to build the matrices and `solve_time` that gives
the time needed to solve the eigenvalue problem. As with `keff`, these
variables contain usable values in expressions that appear after the
first zone definition. Times are given in seconds, and refer to the
actual time elapsed and not to the CPU time used. They both may depend
on how busy the operating system scheduler is, whether memory access is
cached or swapped, etc.


Information about spatial distributions are given as functions. The flux
distribution for group $g$ is given in a function named


    flux_g


where $g$ is supposed to be replaced by the actual group number. The
first---usually the fast---group is number one. The functions take as
many arguments as the spatial dimension of the problem. They are
point-wise defined, with the independent variable being the location of
the cell center and the independent variable being the flux of the
corresponding group. Nevertheless, the functions are linearly
interpolated, and as such they can be evaluated---or integrated,
differentiated, etc---at any point.

These functions can be used wherever a function can appear in any valid
expression handled by wasora. For example, if used within a
`PRINT_FUNCTION` there is no need to provide the selected range as by
default a point-wise defined function will be printed just in the
definition points. However, if a range is specified, the function will
be interpolated and evaluated at the positions asked for. Consider for
example a one-group one-dimensional problem with 20 nodes. Then

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    FILE nodes   nodes.dat
    FILE interp  interpolation.dat
    PRINT_FUNCTION FILE nodes   flux_1  
    PRINT_FUNCTION FILE interp  flux_1 MIN 0 MAX 100 STEP 0.1

creates a file `nodes.dat` containing twenty lines with the cell
center-flux pairs, while the other file, `interpolation.dat`, contains
one thousand lines with the intermediate interpolated values in a
ready-to-plot text format. If no power setpoint is given
(section [3.7](#sec:xenon-settings){reference-type="ref"
reference="sec:xenon-settings"}), fluxes are in dimensionless units such
that the mean value---as discussed in
section [2.4.2](#sec:power){reference-type="ref"
reference="sec:power"}---is equal to one. Otherwise, fluxes have units
of inverse squared length and inverse time.


Even though the functions are point-wise defined, they can be evaluated
at any point in space. For example, if there is a neutron detector
located at $\ensuremath\mathbf{r}=(153, 241, 308)$ then its reading may
be assigned to a variable by

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    detector1 = flux_2(153, 241, 308)

although there may not be a cell located at that point. Details about
function interpolation can be found in wasora documentation.


Information about the resulting cell cross sections values can be
obtained by using the functions


    D_g
    SigmaT_g
    SigmaA_g
    nuSigmaF_g
    ESigmaF_g


They are also linear-interpolated point-wise defined functions
containing the macroscopic cross section associated to the cell
discussed in section [2.2](#sec:xs-association){reference-type="ref"
reference="sec:xs-association"}, according to the selected method as
explained in section [3.8](#sec:scheme-settings){reference-type="ref"
reference="sec:scheme-settings"}.


There is one function named


    bc


that is also point-wise defined that shows which cells contain boundary
condition equations instead of the diffusion equation. It is non-zero
only for the location of the cells that have boundary conditions
equations. The actual value is not important, as is an internal
identification number for each type of condition.


There are five additional functions which give information about the
actual original continuous spatial cross sections distributions, i.e.
before the association to the discrete cells:


    cD
    cSigmaT
    cSigmaA
    cnuSigmaF
    cESigmaF


The 'c' in the name means "continuous", and it is included because it is
not so common to ask milonga for the continuous cross section
distribution as having variables called `D` or `nuSigmaF`. Thus, the
extra 'c' avoids name conflicts. Note that these functions do not
incorporate the group number in its name. Rather, they take $N+1$
arguments, being the first one the group number and the rest the spatial
position. For example, to evaluate the total cross section at point
$\ensuremath\mathbf{r}=(10.5,25.2,4.3)$ of group number 3, ask for

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PRINT cSigmaT(3,10.5,25.2,4.3)

They are not point-wise defined, so when using them within a
`PRINT_FUNCTION` keyword, a range is mandatory, as explained below.


In general, one or more of these functions are part the desired output
of milonga. This information might be passed forward to other codes for
further calculation or might be written to files for analysis and
interpretation (as in the last step of
figure [\[fig:uncertainties\]](#fig:uncertainties){reference-type="ref"
reference="fig:uncertainties"}). These features are part of the common
wasora framework, and are explained in the associated documentation.
However, for completeness, a brief summary of output keywords is
presented. In order to write information into a file, first the `FILE`
primary keyword has to used


    FILE id path [STEP]


The two mandatory arguments are an internal identification string and
the actual file path, which can be either absolute or relative to the
actual directory. If the optional secondary keyword `STEP` is given,
then one different file for each step is created, appending the step
number to the file path. In any case, if the file already exists, it is
overwritten.

The keywords used to print data are


    PRINT [FILE id] [CONDITION expr] [expr1 | TEXT string1] [expr2 | TEXT string2] [...]
    PRINT_FUNCTION [FILE id] [CONDITION expr] function1 [function2] [function3] [...] [MIN x_min [y_min [z_min]] MAX x_max [y_max [z_max]] STEP x_step [y_step [z_step]]]


If no file is given, default is to write to standard output. If a
condition is given, the argument expression is evaluated and the print
instruction is executed only if it is nonzero. For the `PRINT` keyword,
a list of expressions or text strings are expected. Each `PRINT` keyword
produces a single line. For the `PRINT_FUNCTION` keyword, a list of
functions is expected. All the requested functions ought to have the
same number of arguments. The `PRINT_FUNCTION` produces one line for
each requested point, with the values of the independent variables in
the first columns and then one column for each requested function with
the value the function evaluates to at that point. If the first function
is point-wise defined, it is not needed to give the range with the
`MIN`, `MAX` and `STEP` secondary keywords, and the output points are
the definition points of the first function, even if the other functions
are not point-wise defined. Anyway, the output points can be defined
giving minimum, maximum and increment values for the independent
variables, according to the number of arguments of the functions. If the
first function is not point-wise defined (for example given by an
algebraic expression) the the range is mandatory. Actually, these
primary keywords can take more optional secondary keywords. See the
wasora documentation for further information.

For example, valid constructions for the `PRINT` keyword are

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PRINT 1+1
    PRINT TEXT "the eigevalue problem insumed" solve_time TEXT "seconds"
    FILE output outputs/results.txt
    PRINT FILE output   (keff-1)/keff  solve_time
    PRINT CONDITION greater(keff,1) TEXT "the reactor is supercritial" 

To produce a text file containing the fast and thermal flux
distribution, simply enter

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    FILE flux flux.dat
    PRINT_FUNCTION FILE flux  flux_1 flux_2

To obtain the continuous macroscopic cross sections, a range is
mandatory because they are not point-wise defined. For example, in a
two-dimensional problem, one may want to do

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PRINT_FUNCTION cSigmaA_1 cnuSigmaF_2 MIN 1 0 0 MAX 1 x_bare_length y_bare_length STEP 1 0.1*x_bare_length/x_cells 0.1*y_bare_length/x_cells


to obtain a representation of the absorption and nu-fission cross
sections with a resolution of $10\times10$ grid for each cell.
Chapter [4](#cap:examples){reference-type="ref"
reference="cap:examples"} has a lot of further examples and applications
of these two keywords.

Power and xenon distribution {#sec:xenon-settings}
----------------------------

By default, milonga normalizes the eigenvector obtained by solving the
discrete multigroup diffusion problem to obtain a dimensionless flux
distribution with mean value equal to one unit of flux, as discussed in
section [2.4.2](#sec:power){reference-type="ref" reference="sec:power"}.
To obtain a dimensional flux distribution with physical meaning, a power
setpoint has to be given. Also, values for $E\Sigma_f$ for at least one
zone has to be given, in order to compute the normalization factor. The
power setpoint is entered by means of the special variable


    power


The units of the setpoint power given in variable `power` depend on the
number of spatial dimensions of the problem. For three-dimensional
cases, units are power, i.e., energy and inverse time. In two
dimensions, `power` is a power density per unit length perpendicular to
the $x$-$y$ plane in energy inverse time inverse length. In one
dimension, it is a density per unit area perpendicular to the $x$ axis
in units of energy inverse time inverse squared length.

If the `power` variable is set to a non-zero value, two functions of
space are defined:


    power_density
    xenon


The first one gives the power density dissipated locally at
position $\ensuremath\mathbf{r}$ in units of energy inverse times
inverse cubed length. As with the fluxes, this function is
linearly-interpolated point-wise defined, and can be used to coupled
milonga with thermalhydraulic or plant codes---as with any other
function or variable within the wasora framework. It can also be used
the power peak factor, for example in one dimension:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    peak = max(power_density(x), x, 0, x_bare_length)/(power/x_bare_length)


The other function is called `xenon`, and contains the steady-state
xenon concentration distribution associated to the last computed neutron
flux distribution, in units of inverse cubed length. To compute this
function, equation [\[eq:xenon\]](#eq:xenon){reference-type="eqref"
reference="eq:xenon"} is used. There are some defaults for the yields,
the xenon decay constant and the number of neutrons emitted in one
fission, but they can be changed by means of the following variables


    gamma_I
    gamma_Xe
    lambda_I
    lambda_Xe
    nu


whose default values are

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    gamma_I = 6.386e-2
    gamma_Xe = 2.28e-3
    lambda_I = 2.878e-5    # [1/sec]
    lambda_Xe = 2.10e-5    # [1/sec]
    nu = 2.4

Note that the steady-state xenon concentration does not depend
on $\lambda_\text{I}$, so this variable is not used in this version.
Iodine and xenon yields---i.e. how many nuclei of $^{135}$I and of
$^{135}$Xe appear promptly by a single fission---and $\nu$ may be
changed according to the type of fuel used in the problem. The decay
constants---whose default value is given in inverse seconds---are
physical constants rather than variables, but the user may want to
change $\lambda_\text{Xe}$ to study and try to understand how xenon
poisoning works while in college or while having to develop a power
distribution control algorithm, or to use a time unit different than
seconds.


The microscopic neutron absorption of $^{135}$Xe is entered as a vector
of size $G$ named


    sigmaAXe


whose individual elements can be set as with the `chi` vector. The
default behavior is to have zero absorption for every group except for
the last one---considered the thermal group. Assuming a two-group
problem, vector `sigmaAXe` is set to

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    sigmaAXe_1 = 0
    sigmaAXe_2 = 2e-18     # 1/cm^2

that corresponds to a value of two million barns when lengths are
measured in centimeters.


If the nuclear parameters given in the `MATERIALS` keywords do not
depend on the xenon distribution, then the function `xenon` is only
informative and does not need further treating. However, if---as
usual---macroscopic cross sections depend on the xenon distribution by
means of direct reference to the `xenon` function, an iterative
calculation is needed.

As discussed in page , xenon and power distributions depend on each
other. In the first calculation, it is assumed that there is no xenon in
the reactor core. If xenon effects are to be taken into account, a new
power calculation is needed, giving rise to a new xenon concentration,
and so on. To tell milonga to make an iterative calculation, the special
variable


    static_iterations


has to be set to a value greater than one. The current iteration number
can be read from the content of the special variable


    static_step


Convergence in the sense of
equation [\[eq:iterations\_convergence\]](#eq:iterations_convergence){reference-type="eqref"
reference="eq:iterations_convergence"} is not checked for. The user has
to estimate a reasonable value for `static_iterations`. However, an
explicit converge check can be made from within the input file, setting
the special wasora variable `done` to non-zero when a certain condition
holds:

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    # [...]
    # problem definition
    # [...]

    static_iterations = 20

    # keep the keff found in the previous step
    kold = keff

    # define zones here so the the new keff is computed at this point
    ZONE fuel MATERIAL fuel

    # give the user some feedback in standard output
    PRINT static_step keff

    # set done to true if convergence is attained
    done = less(abs(keff-kold), 1e-6)


When printing distributions in an iterative calculation, output files
should be defined using the secondary keyword `STEP` to the primary
keyword `FILE`. See the examples in the next chapter.

Numerical scheme settings {#sec:scheme-settings}
-------------------------

The basic numerical scheme and selections about how to associate cross
sections to cells and how to evaluate the leakage term are given by
using the primary keyword `SCHEME` along with its secondary keywords:

\lstset{language=mil,backgroundcolor=\color{white},frame=none}
    SCHEME [VOLUMES | DIFFERENCES] [XS_MEAN | XS_CENTER] [[D_EPSILON | D_MEAN] | [GRAD_D_LOCAL | GRAD_D_NEIGHBORS]] [NO_NEGATIVE_FLUX | ALLOW_NEGATIVE_FLUX]

If the `VOLUMES` keyword is given, the leakage term is treated with the
finite-volumes based method as described in
section [2.1.2](#sec:volumes){reference-type="ref"
reference="sec:volumes"}, while `DIFFERENCES` keywords selects the
finite-differences method of
section [\[sec:differences\]](#sec:differences){reference-type="eqref"
reference="sec:differences"}. `XS_MEAN` associates cross sections to
cells according to
equation [\[eq:cell-xs-mean\]](#eq:cell-xs-mean){reference-type="eqref"
reference="eq:cell-xs-mean"} and `XS_CENTER` uses
equation [\[eq:cell-xs-center\]](#eq:cell-xs-center){reference-type="eqref"
reference="eq:cell-xs-center"}. For finite volumes, `D_EPSILON` computes
the diffusion coefficients at a distance $\epsilon$ of the cell border,
i.e. $D(x\pm\epsilon)$, while `D_MEAN` uses mean values in the proper
half of the cell, i.e. $\left\langle D(x\pm\epsilon)\right\rangle$. For
finite differences, `GRAD_D_LOCAL` makes milonga to evaluate the
gradient of the diffusion coefficient in
equation [\[eq:diff-grad\]](#eq:diff-grad){reference-type="eqref"
reference="eq:diff-grad"} by locally differentiating the continuous
expressions of $D(\ensuremath\mathbf{r})$. On the other hand,
`GRAD_D_NEIGHBORS` computes the gradient by a finite-difference
approximation using first neighbors. When keyword `NO_NEGATIVE_FLUX` is
given, if a flux happens to be negative because of the application of
discrete boundary conditions
(figure [\[fig:bc-planarneg\]](#fig:bc-planarneg){reference-type="ref"
reference="fig:bc-planarneg"}), it is forced to be zero. Otherwise, if
keyword `ALLOW_NEGATIVE_FLUX` is entered, negative fluxes are kept as
computed by the eigenvalue problem solver.

If no `SCHEME` keyword is entered, the default behavior is as if

\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    SCHEME VOLUMES XS_MEAN D_MEAN NO_NEGATIVE_FLUX

was entered. For finite differences, default is `GRAD_D_NEIGHBORS`.


The value of $\epsilon$ is taken from the special variable

    epsilon

that has units of length and is set to $10^{-2}$ by default.

When using `XS_MEAN`, the integral in the numerator of
equation [\[eq:cell-xs-mean\]](#eq:cell-xs-mean){reference-type="eqref"
reference="eq:cell-xs-mean"} is computed by an adaptive method, whose
convergence criteria can be modified by setting the special variable

    xs_mean_tolerance

that is the relative error allowed. The default value is $10^{-2}$
meaning 1% tolerance. The computational effort needed to build the
problem matrices depends heavily on this parameter.


Another keyword that controls the behavior of the solution is

\lstset{language=mil,backgroundcolor=\color{white},frame=none}
    SOLVER [METHOD [SLEPC | NONE | user]] [EIGENVALUE_K | EIGENVALUE_INVERSE_K] [TOLERANCE tol]

Method can be either `SLEPC` (default), `NONE` to skip the eigenvalue
problem solution (useful when having problems with the nodalization and
preparing a new problematic input file) or a user-defined routine, as
explained in section [3.14](#sec:user-solver){reference-type="ref"
reference="sec:user-solver"}. If `EIGENVALUE_K` is entered, the
eigenvalue problem defined by
equation [\[eq:direct\_keff\]](#eq:direct_keff){reference-type="eqref"
reference="eq:direct_keff"} is solved (default), whilst
problem [\[eq:inverse\_keff\]](#eq:inverse_keff){reference-type="eqref"
reference="eq:inverse_keff"} is solved for `EIGENVALUE_INVERSE_K`.
Default is "direct $k_\text{eff}$ formulation". The optional tolerance
is the value the relative residual $\| r_\text{rel} \|$ has to take in
order to considered the eigenvalue solution as converged, according to
equation [\[eq:eigen-converged\]](#eq:eigen-converged){reference-type="eqref"
reference="eq:eigen-converged"}. Default is to take SLEPc's default
value that is equal to $10^{-7}$.

Parametric calculations
-----------------------

The capability to perform parametric calculations is a feature of the
wasora framework, and thus just a brief summary is given here. By a
parametric calculation, it is understood a series of computation of the
problem solution with a systematic variation of one or more parameters.
In the case of milonga, one or more variables can be swept through an
interval, either linearly or exponentially.

To define a linear parametric study, the following keyword has to be
entered


    PARAMETRIC variable minimum maximum step


The values `min`, `max` are the ends of the interval to be studied and
`step` is the increment. They cannot be expressions, they have to be
floating point constants. If there selected variable has one or more
explicit assignments in the input file, they are ignored. For example,


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PARAMETRIC a   1  10  1
    a = 5
    PRINT a


will output the natural numbers up to ten instead of ten values of five.

If more than one `PARAMETRIC` keyword is given, multidimensional studies
can be performed. The variable in the first `PARAMETRIC` keyword is
incremented up to the maximum value, then the variable in the second
`PARAMETRIC` keyword is incremented one step and the first variable
starts from the minimum value again, and so on.

Instead of incrementing the variable with a fixed value, a
logarithmic---or exponential, depending on the point of view---study can
be done with


    PARAMETRIC_LOG variable minimum maximum step


In this case, `step` is a multiplicative constant by which is multiplied
the selected variable in each calculation.


Each calculation step is independent from the previous one, i.e., the
initial conditions of step $n+1$ are not the last conditions of
step $n$. By default, milonga waits for the previous step to finish
before launching the next one. However, because of this reason, this
kind of calculations are especially suitable to take advantage of
multicore architectures by running simultaneously two or more parametric
steps. To have milonga run in parallel, use


    MAX_DAUGHTERS  n


The parameter `n` is an integer representing the maximum number of
process that may be launched in parallel. If it is equal to zero, the
number of available cores is taken. Default is one, meaning no
parallelization.

When performing parametric calculations, output files are renamed and
the number of step is attached to it. Standard output can be shared by
all the steps, but because of random uncertainties in the operating
system scheduler the order of printing may not be the expected. For
example,


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    MAX_DAUGHTERS  8
    PARAMETRIC a   0  10  1
    PARAMETRIC b   0  10  1
    PRINT a b a+b


may not output the $[1,1]\times[10,10]$ square in the correct order and
plotting utilities may complain about it. To overcome this difficulty,
one solution is to output each line to a file and the concatenate all of
them into a single one in the correct order. First,


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    MAX_DAUGHTERS  8
    PARAMETRIC a   0  10  1
    PARAMETRIC b   0  10  1
    FILE out tmp
    PRINT FILE out a b a+b


and then in the terminal


\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
     $ milonga square.mil
     $ cat tmp.* > output
     $ rm tmp.*
     $

Note that although some calculations may be performed either by using
parametric or iterative calculations, they are essentially different in
their nature. Each instance of a parametric calculation is independent
from previous one---indeed, they are different processes and can be run
in parallel---and all the variables, functions and parameters are
initialized for each instance. On the other hand, each step of a
parametric calculation keeps the variables and functions computed in the
previous step. Therefore, they are not independent and they cannot be
run in parallel.

Debugging and benchmarking
--------------------------

One of the most important issues of an engineering computer code is the
the running time. There are many parameters that affect the
computational time, and more often than not they depend on the
particular problem being solved. Milonga can dump some information about
the numerical scheme that may help the user to reduce the running time
or simply provide extra output.

To obtain extra debug output, use the following keyword


    DEBUG filename [MATRICES_ASCII] [MATRICES_BINARY]


The `filename` argument should be a full or relative path to the text
file with the extra information. Variables `build_time` and `solve_time`
are included in the debug file. If the optional secondary keyword
`MATRICES_ASCII` is given, an ASCII representation of the matrices $R$
and $F$ are included in the file. If the secondary keyword
`MATRICES_BINARY` is given, two binary files called `R` and `F` in the
execution directory are created. They contain a binary dump of the
matrices as created by PETSc's `PetscViewerBinaryOpen` and `MatView`
functions, that can be accessed by using `MatLoad` [@petsc-user-ref].


There is a lot of useful information about the problem being solved that
should be dumped in the debug file. Future versions of milonga will
include further debugging and bechmkarking options.

Commandline arguments replacement
---------------------------------

Milonga takes as a first argument the name of the input file to process.
Using a wasora feature, more optional arguments can be given, and then
replaced in the input. If a the token `$1` is entered in the input file,
it will be replaced by the first non-file commandline argument, `$2`
with the next one, and so on. For example


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PRINT $1


will print whatever argument is given in the commandline after the input
file:

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
     $ milonga calc.mil 1+1
    2.000000e+00
     $ 

This can be useful to provide filenames from the commandline


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    FILE output        $1
    FILE distribution  $2


or to choose between several options without having to modify the input


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    SCHEME $1


If less arguments than the required are provided, milonga fails with an
error.

Advanced eigenvalue solver settings
-----------------------------------

Besides commandline argument replacement, milonga can handle generic
PETSc and SLEPc arguments that are to be passed to their respective
options databases. Thus, run-time options can be changed by using the
standard commandline key mechanism provided by these
libraries [@petsc-user-ref; @slepc-users-manual].

For example, the eigenvalue solver and the preconditioner can be changed

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
     $ milonga test.mil -eps_type power -st_pc_type jacobi

Note that some solvers (and preconditioners) may work only in either the
direct or inverse-$k_\text{eff}$ problem
(section [3.8](#sec:scheme-settings){reference-type="ref"
reference="sec:scheme-settings"}). See the PETSc and SLEPc
manuals [@petsc-user-ref; @slepc-users-manual] for further reference.

External coupling
-----------------

To be explained.

User-defined eigenvalue solver {#sec:user-solver}
------------------------------

To be explained.

If anything goes wrong
----------------------

As stated in section [1.2](#sec:basis){reference-type="ref"
reference="sec:basis"}, one of milonga's the design basis vectors is
"resisting" inconsistent or non-physical values of cross sections. This
may lead to singular matrices, zero pivots or other mathematical
problems involving invalid operations such as divisions by zero.
Depending on the particular situation, one may want the code to behave
differently. For example, sometimes when a singularity is found it would
be desirable to tell the user that there are some problems in the
physics entered in the input and to stop the calculation. But when
running parametrically or inside an optimization loop, it would be
convenient just to ignore those cases and proceed with the next
iterations.

All mathematical operations performed by wasora that may give rise to an
invalid result---division, square roots, logarithms, etc---are checked
and, in case the result is not a number, a handler is called. What the
code does is controlled by the following optional keyword in the input
file


    ON_NAN [INFORM | BE_QUIET] [QUIT | CONTINUE]


Default behavior is to print a message to the standard error output
(`INFORM`) and quit immediately (`QUIT`), but any combination of giving
the error or not and continuing the execution or quitting may be given.

As wasora uses some features from the GNU Scientific
Library [@gsl-manual], there may be run-time errors when accessing one
of these routines---for example divergence of an integral, inconsistent
data for interpolation, etc. When an error is raised during a GSL call,
also the default behavior is to report it to the standard error and
quit. This can be changed using

    ON_GSL_ERROR [INFORM | BE_QUIET] [QUIT | CONTINUE]

Finally, when working with milonga, a PETSc error may arise if the
neutronic problem is not well defined. Again, the default behavior is to
report and quit. The keyword to modify this is

    ON_ERROR [INFORM | BE_QUIET] [QUIT | CONTINUE]

For example,


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    PRINT 1/0


gives


\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    error: NaN found


while


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    ON_NAN CONTINUE BE_QUIET
    PRINT 1/0


results in

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    inf


In the same sense, for example if a zone with a diffusion coefficient
equal to zero is entered, probably PETSc will complain through milonga
with something like


\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    error: PETSc error 71-1 'Zero pivot row 99 value 3.20577e-15 tolerance 1e-12' in src/mat/impls/aij/seq/aijfact.c MatLUFactorNumeric_SeqAIJ:574


Note that the keywords `ON_` define the overall behavior. They do not
work in a by-block basis. The last keyword entered is how the code will
behave for the whole input. For example


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    ON_NAN CONTINUE BE_QUIET
    PRINT 1/0
    ON_NAN INFORM
    PRINT sqrt(-1)


gives


\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    error: NaN found
    inf
    error: NaN found
    0.000000e+00


i.e. as if


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    ON_NAN CONTINUE INFORM
    PRINT 1/0
    PRINT sqrt(-1)


was entered. A final example illustrating `ON_GSL_ERROR`:


\lstset{language=mil,backgroundcolor=\color{mil_fondo},frame=single}
    ON_GSL_ERROR INFORM   QUIT
    PRINT integral(1/(x-1),x,0,1.5)



\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    error: gsl error 21 'bad integrand behavior found in the integration interval' in qag.c

\bibliographystyle{unsrt}
Examples {#cap:examples}
========


\hspace{\fill}
\sf 
\footnotesize
--- ?'Que hacés, que hacés?\
Che, !'Rafael! ?'A dónde vas?\
--- Al cabaret, a milonguear...\
ando buscando una mina\
!'y no la puedo encontrar!

*Retintín (milonga), 1918*


Not only do these examples show what kind of problems can be tackled by
milonga, but also---and more important---how to solve them. To obtain
the most from this chapter, a thorough comprehension of
chapters [2](#cap:equations){reference-type="ref"
reference="cap:equations"} and [3](#cap:input){reference-type="ref"
reference="cap:input"} is recommended. However, as artificial neural
networks do, people learn mainly by examples so probably this chapter
can be read as it was a standalone document.

Examples are grouped by kinds of problems and, with some luck,
introduced in an increasingly complexity order. Even though they are
based on basic reactor physics ideas, they do not refer nor apply to any
real engineering condition. It should be stressed again that these
examples provide the results obtained by solving the neutron diffusion
equation as discussed in
chapter [\[cap:equations\]](#cap:equations){reference-type="eqref"
reference="cap:equations"} and by no means represent to any real nuclear
reactor situation.

In each case, special emphasis is given to the feature that the example
is trying to illustrate, trying to keep as simple as possible the rest
of the problem conditions. Thus, most of the cases are
one-dimensional---except those that illustrate how to manage higher
dimensions problems---in order to avoid introducing extra complications.


The output information and figures that conform the result of the
problems solved in this chapter were generated automatically from a
series of scripts. Therefore, all the examples have the same structure.
First, a brief description of the problem along with a figure is given.
Then, there may be one or more milonga runs for each case, each having
an explanatory paragraph. Then, the actual input file is shown, along
with a view of a fictitious terminal where the invocation of milonga and
some further processing for generating graphical plots and results is
given. In fact, the actual figures shown in this chapter were obtained
by executing the commands shown in the terminal by using the
free-as-in-beer-but-not-as-in-speech plotting program
gnuplot [@gnuplot]. All figures are vector-based, so they can be zoomed
in without loosing resolution. A milonga executable with debugging
symbols was used to run the examples shown, so running times may differ
when using optimized versions. The milonga and gnuplot input files are
available in the distribution package, so there is no need to cut &
paste from this document to further study, modify and tweak the
examples. The commands shown in the sample terminal outputs should be
reproducible using the provided files. Feedback about these examples or
other desired test cases are welcome.



Enough chitchat. It is time to *milonguear*!

Cases with analytical solutions
-------------------------------

The examples shown in this first section have two advantages that are
somehow related to each other: they are both fairly simple and have an
analytical solution. Thus, they are suitable for introducing the usage
of milonga and at the same time can be used to test whether things are
being done well or not, and moreover, to what extent.

Most of the them are one-group cases, as few-group problems either do
not have analytical solutions or they are so complicated that do not
serve the point of being a reference. The only two-group problem shown
here corresponds to an infinite reactor. The rest of the problems are
bare one-speed problems, except for the last one that represents a
two-zone slab reactor.

### Homogeneous bare slab

This problem consists of an homogeneous one-dimensional bare reactor
with constant one-group neutron cross sections, whose solution by the
diffusion approximation is one of the first examples of the solution of
the diffusion equation in nuclear reactor theory courses (see for
example [@lamarsh section 9-2] and [@duderstadt section 5.III.C]). In
particular a one-meter width slab is considered, with the one-group
constants shown in the figure.

![image](examples//01-analytical/01-bare_slab/bareslab.pdf)

#### keff.mil {#keff.mil .unnumbered .unnumbered}

This example computes the the effective multiplication factor
$k_\text{eff}$ for the slab above. See the comments in the input file
for details about the implementation in milonga. As can be seen, a small
number of keywords are used and thus most of the default options are
automatically selected. The only output obtain in the terminal after the
execution of the code is the multiplication factor as computed from the
numerical eigenvalue problem in the default scientific notation format.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # computation of the keff for a one-dimensional
    # one-group bare slab with constant XS

    # define the number of spatial dimensions and groups
    PROBLEM DIMENSIONS 1 GROUPS 1

    x_bare_length = 100           # bare length in cm
    x_cells = 100                 # number of cells

    ZONE fuel MATERIAL fuel  # a single zone spanning the whole bare length
                             # composed of material fuel, defined below

    # definition of the fuel XS for the one-group units are those of
    # the inverse x_bare_length, i.e. 1/cm
    MATERIAL fuel {
      D_1        1
      SigmaA_1   0.010
      nuSigmaF_1 0.012
    }

    # print the multiplication factor in the standard output
    # with the default format, i.e. scientific notation
    PRINT keff

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga keff.mil
    1.092211e+00	
    $

#### comparisson.mil {#comparisson.mil .unnumbered .unnumbered}

Taking advantage milonga's---actually wasora's---mathematical
capabilities, it is fairly easy to compare the numerical result obtained
by using a spatial discretization with 100 cells to the actual
analytical solution for the one-dimensional one-group bare reactor. In
this case, the input for the numerical problem is almost the same, but
the cross sections are entered as variables and used both for the
definition of the material and for the computation of the analytical
value of the effective multiplication factor, which is

$$k_\text{eff} = \frac{\nu\Sigma_f}{\displaystyle \Sigma_a + D \left(\frac{\pi}{a}\right)^2}$$
where $a$ is the slab width. The output is now extended to include an
informative string and a decimal format is used to show at which digit
the two solutions start to differ. The absolute difference is shown in
scientific notation. In the comments, an explanation is given about the
order in which the algebraic expressions are computed. If you are
tempted to change the number of cells to see what happens with the
difference, please do it now. Nevertheless, in the examples of
parametric calculations you will get plenty of the kind of fun you are
searching for.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}

    # comparisson between the numerical and analytical solutions
    # for the one-dimensional bare slab

    PROBLEM DIMENSIONS 1 GROUPS 1
    x_bare_length = 100
    x_cells = 100

    # these are plain variables used both for entering the XS in
    # the MATERIAL keyword and in computing the analytical solution
    D = 1
    SigmaA = 0.01
    nuSigmaF = 0.012

    # the place where the MATERIAL keywords appear in the input
    # is irrelevant, but the place where the ZONE keywords appear
    # is not. All the algebraic expressions that appear before the
    # first ZONE are computed before solving the eigenvalue problem,
    # and the rest are computed after solving the eigenvalue problem.
    # This is important in this example because the XS are evaluated
    # with the information available before the first ZONE keyword.
    # Thus, the variables containing the XS have to be evaluated
    # before the following line.
    ZONE fuel MATERIAL fuel

    # the XS can be any valid mathematical expression, including
    # constants, variables, algebraic expressions, functions, etc.
    # Further examples of XS will be given in subsequent examples.
    # In this one, we are using the variables defined above, that
    # should be "filled" before the first ZONE keyword
    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaF_1 nuSigmaF
    }


    # the analytical solution can be evaluated now or before the
    # ZONE, it does not matter. But, if the keff variable is
    # refered to before the ZONE, its value is equal to one.
    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + D*(pi/x_bare_length)^2)


    # print the solutions
    PRINT TEXT "analytical keff = " k    FORMAT %.8lf
    PRINT TEXT "numerical  keff = " keff FORMAT %.8lf
    PRINT TEXT "     difference = " abs(k-keff)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga comparisson.mil
    analytical keff = 	1.09220381	
    numerical  keff = 	1.09221091	
         difference = 	7.101524e-06	
    $

#### flux.mil {#flux.mil .unnumbered .unnumbered}

The last example regarding the bare slab shows how to easily plot a
function of a single variable. Even though the flux is computed at the
cells centers, milonga provides a continuous function by interpolating
the results that can be integrated, differentiated, etc. This example
shows how to plot the flux distribution at the computed points using the
`PRINT_FUNCTION` without the necessity of giving the actual plot points.
Other examples will deal with the continuous interpolation.

This example also shows how to compare the numerical solution with the
analytical one (remember the obtained flux has mean value equal to one):

$$\phi(x) = \frac{\pi}{2} \sin\left( \frac{\pi}{a} \cdot x \right)$$ for
one hundred cells (examples showing how this difference varies with the
cell number are about to come). Milonga gives as output three columns
containing the $x$-coordinate, the numerical flux and the analytical
flux. This information can be easily plotted using gnuplot, as shown in
the terminal output by catting the files `fluxes.gnuplot` and
`diff.gnuplot`. These are listed here for completeness. In other
examples the files used to obtain the output figures may or may not be
shown depending on the complexity of the figure. However, the files
generated by milonga can be plotted or further processed by any other
software as they consist of tab-separated columns. In this example the
flux distribution is written to the standard output, and is redirected
to a file from the commandline as shown in the terminal output.

Note that in this example, the two plotted columns do not refer to the
same quantities, as the first is the mean flux in each cell and the
former is the continuous flux evaluated at the cell center, as discussed
in the mathematical development of the finite volumes scheme. However,
the difference should be minimal and they are assumed to be equal as in
equation [\[eq:approx-meanflux-centerflux\]](#eq:approx-meanflux-centerflux){reference-type="eqref"
reference="eq:approx-meanflux-centerflux"}.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # one-dimensional bare slab flux distribution

    PROBLEM DIMENSIONS 1 GROUPS 1
    x_bare_length = 100
    x_cells = 100
    D = 1
    SigmaA = 0.01
    nuSigmaF = 0.012

    ZONE fuel MATERIAL fuel

    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaF_1 nuSigmaF
    }

    # the flux is given in a function of the variable x called
    # flux_1(x). It is defined for every value of x by a linear
    # interpolation of the fluxes computed at each cell. By
    # default, the flux is such that the mean value is equal to
    # one. On the other hand, an algebraic function with the
    # analytical solution of the problem  having a mean value
    # equal to one can be constructed as follows:
    pi = 4*atan(1)
    FUNCTION phi(x) = pi/2*sin(pi*x/x_bare_length)

    # as flux_1 is defined in a finite number of abscissae,
    # by default it is printed in these values of x if no
    # range is given. Function phi is algebraic, but it uses
    # the range of the first function in the list.
    PRINT_FUNCTION flux_1 phi

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga flux.mil > flux.dat
    $ cat fluxes.gnuplot
    set terminal pdf
    set output "fluxes.pdf"
    set title "bare-slab flux"
    set xlabel "x"
    set border 3
    set tics nomirror
    set tics scale 0.2
    set ytics format "%.1e"
    plot "flux.dat" lt 3 ps 0.5 ti "numerical", "flux.dat" u 1:3 lt 7 w l ti "analytical"
    $ gnuplot fluxes.gnuplot
    $ cat diff.gnuplot
    set terminal pdf
    set output "diff.pdf"
    set title "relative error between theoretical and numerical fluxes"
    set xlabel "x"
    set border 3
    set tics nomirror
    set tics scale 0.2
    set ytics format "%.1e"
    set ytics 1e-4
    unset key
    plot "flux.dat" u 1:(($3-$2)/$3) w lp pt 2 ps 0.25
    $ gnuplot diff.gnuplot
    $

![image](examples//01-analytical/01-bare_slab/fluxes.pdf)

![image](examples//01-analytical/01-bare_slab/diff.pdf)

### Homogenoeus bare square

This example is an extension to two dimensions to the bare slab. It
consists of an homogeneous one-speed bare square. The gnuplot input
files are also shown in the terminal output to illustrate how to obtain
two-dimensional contour levels for the flux. The first example has the
same number of cells in each direction, while the second one does not.
Also, an illustration of how to use both the standard output and output
files using the `FILE` keyword is given.

![image](examples//01-analytical/02-bare_square/baresquare.pdf)

#### uniform.mil {#uniform.mil .unnumbered .unnumbered}

The input that follows solves the one-group neutron diffusion equation
in an homogeneous square, where the variable `y_bare_length` is set
equal to `x_bare_length`. The number of cells in each dimension is the
same, by also setting the variable `y_cells` equal to `x_cells`. The
$k_\text{eff}$ is evaluated both analytically and numerically, and the
flux distribution is written into a file that can be directly plotted
with gnuplot or another similar tool.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # two-dimensional bare slab flux distribution

    PROBLEM DIMENSIONS 2 GROUPS 1
    x_bare_length = 100
    y_bare_length = x_bare_length  # this forces the geometry to be a square
    x_cells = 100
    y_cells = x_cells              # this forces cells to be squares

    D = 1
    SigmaA = 0.01
    nuSigmaF = 0.012

    ZONE fuel MATERIAL fuel

    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaF_1 nuSigmaF
    }

    # analytical effective multiplication factor
    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + D*((pi/x_bare_length)^2 + (pi/y_bare_length)^2))

    # write the flux distribution to a file
    FILE flux uniform.dat
    # the flux is evaluated at the computed points automatically
    PRINT_FUNCTION FILE flux flux_1

    # and the keff to the screen
    PRINT TEXT "  numerical keff =" keff
    PRINT TEXT " analytical keff =" k
    PRINT TEXT "      difference =" abs(keff-k)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga uniform.mil
      numerical keff =	1.002190e+00	
     analytical keff =	1.002178e+00	
          difference =	1.195821e-05	
    $ cat uniform2d.gnuplot
    set terminal pdf
    set output "uniform2d.pdf"
    set title "flux with an uniform mesh"
    set size square
    plot "uniform.dat" w image
    $ gnuplot uniform2d.gnuplot
    $ cat uniform3d.gnuplot
    set terminal pdf
    set output "uniform3d.pdf"
    set title "flux with an uniform mesh"
    set ticslevel 0
    unset key
    splot "uniform.dat" u 1:2:3 palette pt 59 ps 0.1
    $ gnuplot uniform3d.gnuplot
    $

![image](examples//01-analytical/02-bare_square/uniform2d.pdf)

![image](examples//01-analytical/02-bare_square/uniform3d.pdf)

#### nonuniform.mil {#nonuniform.mil .unnumbered .unnumbered}

This example is quite similar to the previous one, but the number of
cells in the $y$-direction is one-fifth of the number of cells in the
$x$-direction. The flux distribution output is evaluated at the computed
points automatically by milonga. The effects of the different
discretizations is studied in other examples.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # one-dimensional bare slab flux distribution

    PROBLEM DIMENSIONS 2 GROUPS 1
    x_bare_length = 100
    y_bare_length = x_bare_length
    x_cells = 100
    # the following line sets the number of y_cells to
    # one fifth of whatever x_cells is set to. If the
    # result is not integral, the value is rounded to
    # the nearest integer
    y_cells = x_cells/5 

    D = 1
    SigmaA = 0.01
    nuSigmaF = 0.012

    ZONE fuel MATERIAL fuel

    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaF_1 nuSigmaF
    }

    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + D*((pi/x_bare_length)^2 + (pi/y_bare_length)^2))

    FILE flux nonuniform.dat
    # no need to change this line with respect to the previous
    # example, milonga knows where function flux_1 is defined
    PRINT_FUNCTION FILE flux flux_1

    PRINT TEXT "  numerical keff =" keff
    PRINT TEXT " analytical keff =" k
    PRINT TEXT "      difference =" abs(keff-k)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga nonuniform.mil
      numerical keff =	1.002253e+00	
     analytical keff =	1.002178e+00	
          difference =	7.522534e-05	
    $ gnuplot nonuniform2d.gnuplot
    $ gnuplot nonuniform3d.gnuplot
    $

![image](examples//01-analytical/02-bare_square/nonuniform2d.pdf)

![image](examples//01-analytical/02-bare_square/nonuniform3d.pdf)

### Homogeneous bare circle

Now a bare circular one-speed reactor is solved. To define the geometry,
secondary keywords have to be given to the primary keyword `ZONE` to
prevent the zone to expand to the whole bare length. Note that the error
committed in the computation of $k_\text{eff}$ is greater than that of
the bare square because of the introduction of discrete circular
boundary conditions. Some capabilities of the wasora's common framework
for working with multidimensional functions are also shown. The
analytical solution involves the use of the lowest-order two-dimensional
Laplacian eigenfunction with null boundary conditions over a circle,
namely Bessel's $J_0$ function. As wasora does not provide evaluation of
Bessel's functions, the profile comparison is made by making use of
gnuplot capabilities.

![image](examples//01-analytical/03-bare_circle/barecircle.pdf)

#### full.mil {#full.mil .unnumbered .unnumbered}

This input computes the one-group flux inside a a full homogeneous
circle of radius $R=100 \, \text{cm}$. The cells that are marked to have
a boundary condition instead of the diffusion equations are shown by
printing the `bc` function. This example uses the default behavior of
forcing negative fluxes to zero at boundary condition cells, effect that
can be seen as a discontinuity in the colors when zooming into the
border of the flux distribution contour plot. The radial flux profile is
computed by defining a one-dimensional function `bessel` from the
two-dimensional flux distribution. It is compared to the analytical
solution by plotting this continuous function along with gnuplot's
Bessel function implementation.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # two-dimensional bare circle

    PROBLEM DIMENSIONS 2 GROUPS 1

    # files to output the full 2D distribution and the radial
    # distribution to compare to the analytical profile. The
    # location of this definition inside the input file does
    # not matter, only the PRINT instructions are executed
    # in the given order
    FILE flux   full.dat
    FILE radial radial.dat

    # square bare geometry and square cells
    x_bare_length = 100
    y_bare_length = x_bare_length
    x_cells = 100
    y_cells = x_cells

    D = 1
    SigmaA = 0.01
    nuSigmaF = 0.012

    # only one zone is defined, a centered circle of radius equal
    # to half the bare length. Cells that lie outside the circle
    # are forced to have zero flux
    ZONE fuel MATERIAL fuel X_CENTER x_bare_length/2 Y_CENTER x_bare_length/2 RADIUS x_bare_length/2

    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaf_1 nuSigmaF
    }

    # the analytical multiplication factor
    nu0 = 2.4048    # first zero of bessel's j0 function
    k = nuSigmaF/(SigmaA + D*(nu0/(x_bare_length/2))^2)

    # definition of a one-dimensional function equal to the flux
    # distribution evaluated in an horizontal radius, i.e. the
    # radial profile to compare to the bessel solution
    FUNCTION bessel(r) = flux_1(r+x_bare_length/2, x_bare_length/2)

    # write the flux distribution and which cells have boundary
    # conditions to one file
    PRINT_FUNCTION FILE flux   flux_1 bc
    # and the radial profile to another one. Note that bessel(r) is
    # an algebraic function, so a range is mandatory. As the step is
    # less than delta_x, this is an interpolation of the actual flux.
    PRINT_FUNCTION FILE radial bessel MIN 0 MAX x_bare_length/2 STEP 0.5*x_bare_length/x_cells

    PRINT TEXT "  numerical keff =" keff
    PRINT TEXT " analytical keff =" k
    PRINT TEXT "      difference =" abs(keff-k)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga full.mil
      numerical keff =	9.773084e-01	
     analytical keff =	9.745619e-01	
          difference =	2.746505e-03	
    $ gnuplot full2d.gnuplot
    $ gnuplot full3d.gnuplot
    $ gnuplot fullbc.gnuplot
    $ cat radial.gnuplot
    set terminal pdf
    set output "radial.pdf"
    set title "radial flux shape vs. analytical solution"
    set xzeroaxis
    set border 3
    set tics scale 0.2
    set tics nomirror
    set yrange [-0.1:3]
    set key bottom left Left reverse
    plot "radial.dat" ps 0.5 lt 3 ti "numerical", 2.914*besj0(x/50*2.4048) lt 7
    $ gnuplot radial.gnuplot
    $

![image](examples//01-analytical/03-bare_circle/fullbc.pdf)

![image](examples//01-analytical/03-bare_circle/full2d.pdf)

![image](examples//01-analytical/03-bare_circle/full3d.pdf)

![image](examples//01-analytical/03-bare_circle/radial.pdf)

### Infinite reactor

An infinite reactor is an imaginary reactor that spans the whole
universe and thus, the neutron flux does not depend on the position. For
the multigroup formulation, there is analytical expression for the
multiplicative factor, called $k_\infty$. This example uses a trick to
fool milonga into thinking it has to solve an two-group infinite reactor
and illustrates how to obtain extra information about the problem.

#### kinf.mil {#kinf.mil .unnumbered .unnumbered}

First, the analytical solution for a two-group infinite reactor with no
upscattering

$$k_\infty = \frac{\nu\Sigma_{f1}}{\Sigma_{a1}} + \frac{\nu\Sigma_{f2}}{\Sigma_{a2}} \cdot \frac{\Sigma_{12}}{\Sigma_{a1}}$$

is computed and stored into a variable. Then, the numerical solution is
computed. It involves a one-dimensional two-group problem with three
spatial cells. The boundary conditions are set to mirror in both ends,
thus forcing all the fluxes to be equal as in the infinite reactor. This
is the first example of a two-group energy problem. See how the
downscattering therm is entered in the `MATERIAL` keyword. The main
output of this example is the textfile with the debugging information,
including an ASCII representation of the matrices of the problem, as
catted in the terminal output.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # computation of the k_infinite value of a medium
    # in a 2-group problem formulation with no upscattering
    #

    # two-group XS as variables
    nuSigmaF1 = 5e-4
    SigmaA1 = 4e-3
    nuSigmaF2 = 5e-3
    SigmaA2 = 9e-3
    SigmaS12 = 8e-3

    # analytical solution
    #########################
    analytical = nuSigmaF1/SigmaA1 + nuSigmaF2/SigmaA2 * SigmaS12/SigmaA1


    # numerical solution
    ##########################
    PROBLEM DIMENSIONS 1 GROUPS 2

    # we need at least three cells: two for the boundary
    # conditions and one for the diffusion equation
    x_cells = 3

    # set both left and right bc to mirror, so all theoretical
    # three cells will end up with the same neutron flux level
    BOUNDARY_CONDITIONS {
     X_MIN MIRROR
     X_MAX MIRROR
    }

    # the bare length does not matter, because of the
    # selected boundary conditions. However, this value must
    # be different from zero
    x_bare_length = 1


    # dump information in a text file, including an ascii
    # representation of the matrices of the problem
    DEBUG infinite.txt MATRICES_ASCII

    ZONE reactor MATERIAL mix

    # the value for the diffusion coeffcient is irrelevant
    # because the laplacian of the flux is indentically zero
    # pay attention to the difference between the keywords
    # and the variables below
    MATERIAL mix {
      D_1        1
      nuSigmaF_1 nuSigmaF1
      SigmaT_1   SigmaA1

      D_2        1
      nuSigmaF_2 nuSigmaF2
      SigmaT_2   SigmaA2

      SigmaS_1.2 SigmaS12
    }

    PRINT TEXT "analytical keff = " analytical            FORMAT %.18f
    PRINT TEXT "numerical keff  = " keff                  FORMAT %.18f
    PRINT TEXT "difference      = " abs(analytical-keff)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga kinf.mil
    analytical keff = 	1.236111111111111160	
    numerical keff  = 	1.236111111110656635	
    difference      = 	4.545253e-13	
    $ cat infinite.txt
    ------------------                 -------------       ------
    milonga 0.1 (linux-i686)
    free nuclear reactor core analysis code
    debugging and benchmarking output
    --------       -------      -----     ----    ---

    execution date is Sat Jul 16 12:58:58 2011


    invocation commandline was

     milonga kinf.mil


    ---------         -------      ------    ----   ---   
      spatial dimensions: 1
                   x_cells: 3
                   delta_x: 3.333333e-01
           energy groups: 2

            problem size: 6


    ----------      ------      -----    ----
    boundary conditions

      x direction
         left: mirror
        right: mirror




    ------------          ----------         ---------        -------     -----    ----
    eigenvalue solution

                                solution method: krylovschur

                number of requested eigenvalues: 1
                number of converged eigenvalues: 6
             number of iterations of the method: 1
      number of linear iterations of the method: 6

      stopping condition: tol=1e-07, maxit=100

            residual norm 
                     |Ax-kBx|_2 = 4.376e-15
            relative error
              |Ax-kBx|_2/|kx|_2 = 1.957e-15
            error estimate
                   |k - k_real| = 1.764e-57


    ----------         ---------         -------      -----    ---
    removal matrix
     +1.00e+00     0     -1.00e+00     0         0         0    
         0     +1.00e+00     0     -1.00e+00     0         0    
     -9.00e+00     0     +1.80e+01     0     -9.00e+00     0    
         0     -9.00e+00 -8.00e-03 +1.80e+01     0     -9.00e+00
         0         0     -1.00e+00     0     +1.00e+00     0    
         0         0         0     -1.00e+00     0     +1.00e+00


    -------------           ----------          ----------        -------
    fission matrix
         0         0         0         0         0         0    
         0         0         0         0         0         0    
         0         0     +5.00e-04 +5.00e-03     0         0    
         0         0         0         0         0         0    
         0         0         0         0         0         0    
         0         0         0         0         0         0    


    -----------          ---------         ---------        --------     ----    
    ------------------                ---------------              --------------            ------
    transcription of input file:
    --------------              -------------           --------       ------     

    # computation of the k_infinite value of a medium
    # in a 2-group problem formulation with no upscattering
    #

    # two-group XS as variables
    nuSigmaF1 = 5e-4
    SigmaA1 = 4e-3
    nuSigmaF2 = 5e-3
    SigmaA2 = 9e-3
    SigmaS12 = 8e-3

    # analytical solution
    #########################
    analytical = nuSigmaF1/SigmaA1 + nuSigmaF2/SigmaA2 * SigmaS12/SigmaA1


    # numerical solution
    ##########################
    PROBLEM DIMENSIONS 1 GROUPS 2

    # we need at least three cells: two for the boundary
    # conditions and one for the diffusion equation
    x_cells = 3

    # set both left and right bc to mirror, so all theoretical
    # three cells will end up with the same neutron flux level
    BOUNDARY_CONDITIONS {
     X_MIN MIRROR
     X_MAX MIRROR
    }

    # the bare length does not matter, because of the
    # selected boundary conditions. However, this value must
    # be different from zero
    x_bare_length = 1


    # dump information in a text file, including an ascii
    # representation of the matrices of the problem
    DEBUG infinite.txt MATRICES_ASCII

    ZONE reactor MATERIAL mix

    # the value for the diffusion coeffcient is irrelevant
    # because the laplacian of the flux is indentically zero
    # pay attention to the difference between the keywords
    # and the variables below
    MATERIAL mix {
      D_1        1
      nuSigmaF_1 nuSigmaF1
      SigmaT_1   SigmaA1

      D_2        1
      nuSigmaF_2 nuSigmaF2
      SigmaT_2   SigmaA2

      SigmaS_1.2 SigmaS12
    }

    PRINT TEXT "analytical keff = " analytical            FORMAT %.18f
    PRINT TEXT "numerical keff  = " keff                  FORMAT %.18f
    PRINT TEXT "difference      = " abs(analytical-keff)
    $

### Slab with continuously-changing properties

The following example consists of a one-speed slab whose nuclear
parameters depend on the axial coordinate as continuous algebraic
functions of $x$ in such a way that the neutron flux can be explicitly
computed also as algebraic functions.

In particular, a non-dimensional slab reactor spanning the range $[0,1]$
subject to null boundary conditions at $x=0$ and mirror boundary
conditions at $x=1$ with the dimensionless cross sections varying as

$$\begin{aligned}
D(x) &= \frac{1}{2} \, x^2 \\
\Sigma_a(x) &= 2 \, \frac{1-x}{2-x} \\
\nu\Sigma_f(x) &= \frac{x^2}{2x - x^2}\end{aligned}$$ gives rise to a
critical reactor whose dimensionless flux spatial distribution---with
mean value equal to one---is

$$\phi(x) = \frac{3}{2} \, (2x - x^2)$$ as can be checked by replacing
the above expressions into

$$\frac{d}{dx} \left( D(x) \cdot \frac{d \phi}{dx} \right) + \left[ \nu\Sigma_f(x) - \Sigma_a(x) \right] \cdot \phi(x) = 0$$

Note that because of stability reasons, the diffusion coefficient has to
be positive. In this case, $D(x) = 0$ for $x=0$, but this value is never
used in the evaluation of the diffusion equation.

#### 20cells.mil {#cells.mil .unnumbered .unnumbered}

This example numerically solves the slab with continuously-changing
properties and compares it to the analytical solution. It illustrates
for the first time how to enter cross-sections as algebraic expressions
of the spatial variables. It also shows how commandline replacement
works. As it may be interesting to solve this problem using different
schemes to compare their behavior, it would be a nice idea to be able to
select options at run-time instead of having them hard-coded in the
input file. The terminal output shows the execution of milonga with four
combinations of scheme options, whose results are redirected to
different files and then plotted into single figures for comparison.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # one-dimensional one-group slab with
    # continously-changing properties
    #
    # the numerical scheme options should be given in
    # the commandline
    #
    # the analytical solution to this problem is
    # phi(x) = 3/2 * ( 2x - x^2 )
    # with keff = 1 (rho = 0)
    #

    PROBLEM DIMENSIONS 1 GROUPS 1

    # to avoid fixing the scheme options in the input file,
    # the options can be given from the commandline
    #
    # $n gets replaced by the n-th commandline argument
    # after the input file. In this case, if less than 3
    # extra arguments are given, milonga quits with an error
    SCHEME $1 $2 $3

    BOUNDARY_CONDITIONS X_MIN NULL X_MAX MIRROR

    x_bare_length = 1   # dimensionless width
    x_cells = 20        # 20 spatial cells

    # small length tending to 0+ to use for the evaluation of
    # the D(x-/+epsilon) coefficients
    epsilon = 1e-2*x_bare_length/x_cells

    ZONE fuel MATERIAL fuel

    # XS given as algebraic expressions of the coordinate x
    MATERIAL fuel {
      D_1         0.5*x^2
      SigmaA_1    2*(1-x)/(2-x)
      nuSigmaF_1  x^2/(2*x-x^2)
    }

    # analytical solution and relative error
    FUNCTION phi(x) = 3/2 * (2*x - x^2)
    FUNCTION error(x) = (flux_1(x) - phi(x))/phi(x)

    # numerical reactivity (should be equal to zero)
    rho = (keff-1)/keff

    # write reactivity to standard output but commented so
    # gnuplot ignores the line. Note that the hash is escaped
    # to avoid milonga ignoring the rest of the line (i.e. it
    # is not a comment for milonga, it is a comment for gnuplot)
    PRINT TEXT "\# rho =" rho

    # write fluxes, errors and XS distributions into stdout
    PRINT_FUNCTION flux_1 phi error D_1 SigmaA_1 nuSigmaF_1

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga 20cells.mil DIFFERENCES XS_CENTER GRAD_D_NEIGHBORS > differences-center-neighbor.dat
    $ milonga 20cells.mil DIFFERENCES XS_MEAN GRAD_D_LOCAL > differences-mean-local.dat
    $ milonga 20cells.mil VOLUMES XS_CENTER D_MEAN > volumes-center-mean.dat
    $ milonga 20cells.mil VOLUMES XS_MEAN D_EPSILON > volumes-mean-epsilon.dat
    $ gnuplot input.gnuplot
    $ gnuplot solution.gnuplot
    $ gnuplot comparisson.gnuplot
    $

![image](examples//01-analytical/05-continuous/input.pdf)

![image](examples//01-analytical/05-continuous/solution.pdf)

![image](examples//01-analytical/05-continuous/comparisson.pdf)

#### 1000cells.mil {#cells.mil-1 .unnumbered .unnumbered}

This example is similar to the previous ones, except that one thousand
cells are used to show that the numerical solution tends to the
analytical one---especially $k_\text{eff}$---as `x_cells`
$\rightarrow \infty$.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}

    PROBLEM DIMENSIONS 1 GROUPS 1

    BOUNDARY_CONDITIONS X_MIN NULL X_MAX MIRROR

    # even more cells can be used, try it!
    x_cells = 1000
    x_bare_length = 1

    ZONE fuel MATERIAL fuel

    MATERIAL fuel {
      D_1         0.5*x^2
      SigmaA_1    2*(1-x)/(2-x)
      nuSigmaF_1  x^2/(2*x-x^2)
    }

    FUNCTION phi(x) = 3/2 * (2*x - x^2)
    FUNCTION error(x) = (flux_1(x) - phi(x))/phi(x)

    rho = (keff-1)/keff

    PRINT TEXT "\# rho =" rho keff
    PRINT_FUNCTION flux_1 error

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga 1000cells.mil > 1000cells.dat
    $ gnuplot 1000cells.gnuplot
    $

![image](examples//01-analytical/05-continuous/1000cells.pdf)

### Two-zone slab {#ex:twozone}

Consider a one-dimensional one-speed two-zone slab of length $b$,
composed of a material $I$ with an infinite multiplication factor less
than one for $0<x<a$ and of a material $II$ with $k_\infty>1$ for
$a<x<b$ (a similar problem is discussed in reference [@cacuci-hebert]).
These two materials have homogeneous properties. It can be shown that
the resulting flux distribution is

$$\phi(x) = 
 \begin{cases}
  \displaystyle \sinh \left[ \sqrt{ \frac{1}{D_\text{I}} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right) } \cdot x \right] & \text{if $0 \le x < a$} \\
\\
  \displaystyle \frac{ \displaystyle \sinh \left[ \sqrt{ \frac{1}{D_\text{I}} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right) } \cdot a \right] }{ \displaystyle \sin \left[ \sqrt{ \frac{1}{D_\text{II}} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}}\right) }  \cdot \left( b - a \right) \right] } \cdot \sin \left[ \sqrt{ \frac{1}{D_\text{II}} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}}\right) }  \cdot \left( b - x \right) \right] & \text{if $b < x \le b$} \\
 \end{cases}$$ where the effective multiplication factor $k_\text{eff}$
is the solution of the critical condition

$$\begin{aligned}
 &\sqrt{ D_\text{I} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right)} \cdot \tan \left[ \sqrt{ \frac{1}{D_\text{II}} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}}\right) }  \cdot \left( b - a \right) \right] + \\
 &\quad \quad \quad \sqrt{ D_\text{II} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}} \right) } \cdot \tanh \left[ \sqrt{ \frac{1}{D_\text{I}} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right) } \cdot a \right] = 0\end{aligned}$$

If the diffusion coefficients are not equal, then there is a
discontinuity in the flux at $x=a$. This problem poses an interesting
benchmark test for nodal methods.

![image](examples//01-analytical/06-two-zone_slab/twozone.pdf)

#### twozone.mil {#twozone.mil .unnumbered .unnumbered}

The two-zone slab problem is solved numerically with $N=200$ cells and
the obtained solution is compared to the analytical flux distribution.
The left half of the slab contains material I and the right half
contains material II. This example shows on the one hand how to include
more than one zone in milonga and on the other some of the mathematical
capabilities that wasora provides in order to obtain the non-trivial
analytical $k_\text{eff}$ and flux distribution $\phi(x)$ as a
continuous algebraic function. The analytical multiplication factor is
computed by using the `root` functional provided by wasora. Because the
critical condition has a rather bad behavior regarding discontinuities
and infinities, the numerical $k_\text{eff}$ is chosen as the initial
guess and a small interval around it is given as the expected solution
range. The analytical solution is written as an ASCII representation of
an algebraic function that can be directly entered into gnuplot. Also, a
text file `flux.dat` is created with both the numerical and the
analytical fluxes evaluated at the cell centers. Note that in this case,
the interface coincides with a cell border. The user is encouraged to
change the number of spatial cells and analyze how the results change.
In the parametric calculations section, this problem is revisited and
how the difference between the numerical and analytical solutions
changes when the interface does not coincide with cell borders.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}

    # this problem consists of a one-dimensional one-speed two-zone slab
    PROBLEM DIMENSIONS 1 GROUPS 1

    # geometry
    ##################################
    a = 50        # width of zone I-start of zone II
    b = 100       # full width

    x_bare_length = b

    x_cells = 100

    # cross sections
    ##################################
    # defined as variables to have them available for
    # computing the analytical keff
    nu1 = 0.010
    nu2 = 0.015
    a1 = 0.015
    a2 = 0.010
    # note the difference between the two diffusion
    # coefficients to stress the discontinuity of the
    # flux gradient
    d1 = 0.5   
    d2 = 1.5

    MATERIAL fuel {
      D_1         d2
      SigmaT_1    a2
      nuSigmaF_1  nu2
    }

    MATERIAL rod {
      D_1         d1
      SigmaT_1    a1
      nuSigmaF_1  nu1
    }

    # zones
    ##################################
    ZONE rod  MATERIAL rod    X_MIN 0   X_MAX a
    ZONE fuel MATERIAL fuel   X_MIN a   X_MAX b

    # analytical multiplication factor
    ##################################
    # interval to look for the analytical k
    eps = 0.05

    # evaluate k from the analytical critical condition
    # taking as an initial guess the numerical keff
    k = keff
    k = root(sqrt(d1*abs(nu1/k-a1)) * tan(sqrt((1/d2)*abs(a2-nu2/k))*(b-a)) + sqrt(d2*abs(a2-nu2/k)) * tanh(sqrt((1/d1)*abs(nu1/k-a1))*a), k, keff-eps, keff+eps)

    # analytical flux distribution                            
    B1 = sqrt(abs(nu1/k-a1)/d1)
    B2 = sqrt(abs(nu2/k-a2)/d2)
    FUNCTION phi1(x) = sinh(B1*x)
    FUNCTION phi2(x) = sinh(B1*a)/sin(B2*(b-a)) * sin(B2*(b-x))
    norm = b/(integral(phi1(x),x,0,a) + integral(phi2(x),x,a,b))
    FUNCTION phi(x) = norm * if(less(x,a), phi1(x), phi2(x))
    FUNCTION error(x) = phi(x)-flux_1(x)

    PRINT TEXT "\# copy and paste the following lines into gnuplot to obtain"
    PRINT TEXT "\# the continuous flux distribution as a function of x"
    PRINT TEXT "numerical_keff  =" keff FORMAT %.8lf
    PRINT TEXT "analytical_keff =" k    FORMAT %.8lf
    PRINT TEXT "phi1(x) = " norm TEXT " * sinh(" B1 TEXT " * x )" SEPARATOR " "
    PRINT TEXT "phi2(x) = " norm TEXT " * sinh(" B1*a TEXT " ) / sin(" B2*(b-a) TEXT ") * sin(" B2 TEXT "*(" b TEXT "-x ) )" SEPARATOR " "
    PRINT TEXT "phi(x) = ( x <" a TEXT ") ? phi1(x) : phi2(x)" SEPARATOR " "

    # print flux_1(x), phi(x) and error(x)
    FILE flux flux.dat
    PRINT_FUNCTION FILE flux flux_1 phi error

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga twozone.mil
    ## copy and paste the following lines into gnuplot to obtain	
    ## the continuous flux distribution as a function of x	
    numerical_keff  =	1.14799085	
    analytical_keff =	1.14779361	
    phi1(x) =  1.316009e-02  * sinh( 1.121395e-01  * x ) 
    phi2(x) =  1.316009e-02  * sinh( 5.606974e+00  ) / sin( 2.261471e+00 ) * sin( 4.522943e-02 *( 1.000000e+02 -x ) ) 
    phi(x) = ( x < 5.000000e+01 ) ? phi1(x) : phi2(x) 
    $ gnuplot twozone.gnuplot
    $

![image](examples//01-analytical/06-two-zone_slab/flux.pdf)

![image](examples//01-analytical/06-two-zone_slab/fluxzoom.pdf)

General problems
----------------

This section presents general cases that do not have analytical solution
but nevertheless should be familiar to the nuclear engineer. They mostly
involve geometries with material boundaries, where the evaluation of
cell cross sections and the leakage term affect the computational effort
needed and the accuracy of the solution in non-trivial ways. No
three-dimensional example is given mainly because the current version of
milonga does not handle correctly cylindrical nor spherical boundary
conditions. Only one and two-dimensional examples are shown.
Nevertheless, most of the multidimensional features of the code can be
illustrated by using two spatial dimensions.

### Reflected slab

After analyzing the bare slab, reactor physics theory courses focus on
the reflected slab. The most important concept to study in these kinds
of problems is the loss of separability between energy and space and the
appearance of the so-called "thermal shoulder" due to the low absorption
and heavy downscattering in the reflector material when using at least
two energy groups. The example that follows solves the classic symmetric
reflected slab shown in the figure, and the next two introduce some
variations to illustrate different aspects.

![image](examples//02-general/01-reflected_slab/reflectedslab.pdf)

#### symmetric.mil {#symmetric.mil .unnumbered .unnumbered}

The reflected slab shown in the figure using two energy groups is
solved. The definition of the materials is given in a separate file and
included in the main input file to shorten its length and to share the
cross sections with other input files. The materials file is shown in
the terminal window.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # symmetric two-group reflected slab

    PROBLEM DIMENSIONS 1 GROUPS 2 

    a = 20     # reflector length
    x_bare_length = 100
    x_cells = 100

    # the XS can be conveniently entered into a separate file
    # and then included in each problem definition input
    INCLUDE materials.mil

    # first, we define a zone that spans the whole slab
    ZONE fuel   MATERIAL fuel
    # and then we add the two reflectors
    # note that the default behavior is to replace the overlapping
    # intervals with the new XS (to have them added use the
    # INCREMENTAL keyword)
    ZONE refl1  MATERIAL reflector X_MIN 0                X_MAX a
    ZONE refl2  MATERIAL reflector X_MIN x_bare_length-a  X_MAX x_bare_length

    # print the two fluxes to the standard output
    PRINT_FUNCTION flux_1 flux_2

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ cat materials.mil
    MATERIAL fuel {
      D_1         1.500
      SigmaT_1    0.010
      nuSigmaF_1  0.000

      D_2         0.400
      SigmaT_2    0.085
      nuSigmaF_2  0.135

      SigmaS_1.2  0.020
    }

    MATERIAL reflector {
      D_1         2.000
      SigmaT_1    0.000

      D_2         0.300
      SigmaT_2    0.010

      SigmaS_1.2  0.040
    }

    MATERIAL rod {
      D_1         1.500
      SigmaT_1    0.010
      nuSigmaF_1  0.000

      D_2         0.400
      SigmaT_2    0.130
      nuSigmaF_2  0.135

      SigmaS_1.2  0.020
    }
    $ milonga symmetric.mil > symmetric.dat
    $ gnuplot symmetric.gnuplot
    $

![image](examples//02-general/01-reflected_slab/symmetric.pdf)

#### asymmetric.mil {#asymmetric.mil .unnumbered .unnumbered}

More often than not, reactors are not symmetric. If they were, there
would be no point in making full-core calculations like the one in the
preceding example. A small variation to the problem above is to study
what happens when the reactor is not symmetrically reflected, where a
full-core calculation is mandatory.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # asymmetric reflected slab

    PROBLEM DIMENSIONS 1 GROUPS 2 

    a = 20     # left reflector length
    b = 10     # right reflector length
    x_bare_length = 100
    x_cells = 100

    INCLUDE materials.mil

    ZONE fuel   MATERIAL fuel
    ZONE refl1  MATERIAL reflector X_MIN 0                X_MAX a
    ZONE refl2  MATERIAL reflector X_MIN x_bare_length-b  X_MAX x_bare_length

    PRINT_FUNCTION flux_1 flux_2

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga asymmetric.mil > asymmetric.dat
    $ gnuplot asymmetric.gnuplot
    $

![image](examples//02-general/01-reflected_slab/asymmetric.pdf)

#### rod.mil {#rod.mil .unnumbered .unnumbered}

A final variation is introduced here where a small interval of the slab
is replaced by a heavy absorber, acting like a control rod. Moreover,
the new material is located in a position whose limits do not coincide
with the cell boundaries. Besides the flux distribution, this example
shows both the total continuous cross section as a function of $x$ and
the cell cross sections. A zoom over the control rod shows what happens
at the cells that contain the material interfaces. As the default
behavior is to use
equation [\[eq:cell-xs-mean\]](#eq:cell-xs-mean){reference-type="eqref"
reference="eq:cell-xs-mean"}, the result is that the cell cross section
is the weighted average of the adjacent materials.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # asymmetric reflected slab with a control rod

    PROBLEM DIMENSIONS 1 GROUPS 2 

    a = 20              # left reflector length
    b = 10              # right reflector length
    rod_x = 40.123      # position of the rod center
    rod_w = 5.456       # total width of the rod

    # note that the cells are 1cm width
    x_bare_length = 100
    x_cells = 100

    INCLUDE materials.mil

    ZONE fuel   MATERIAL fuel
    ZONE refl1  MATERIAL reflector X_MIN 0                X_MAX a
    ZONE refl2  MATERIAL reflector X_MIN x_bare_length-b  X_MAX x_bare_length
    # the control rod_position (in this case the XS are absolute,
    # had they been incremental, the KEYWORD incremental should have
    # been used)
    ZONE rod    MATERIAL rod       X_MIN rod_x-0.5*rod_w X_MAX rod_x+0.5*rod_w

    PRINT_FUNCTION flux_1 flux_2

    # print the continuous thermal SigmaA in contxs.dat
    # and the cell thermal SigmaA in cellxs.dat
    FILE continuous contxs.dat
    FILE cell       cellxs.dat

    # function cSigmaT needs a range and step, the first argument
    # is the energy group and the second one is the x coordinate
    PRINT_FUNCTION cSigmaT MIN 2 0 MAX 2 x_bare_length STEP 1 0.01 FILE continuous

    # function SigmaT_2 does not need an explicit range
    PRINT_FUNCTION SigmaT_2 FILE cell

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga rod.mil > rod.dat
    $ gnuplot rod.gnuplot
    $ gnuplot xs.gnuplot
    $

![image](examples//02-general/01-reflected_slab/rod.pdf)

![image](examples//02-general/01-reflected_slab/xs.pdf)

![image](examples//02-general/01-reflected_slab/xs-zoom.pdf)

### Reflected circle

The next case to study is the reflected circle. The effects that appear
here are pretty much the same than the ones the appeared in the
reflected slab. However, the fact of having a circular geometry embedded
into a rectangular grid gives rise to interesting forms of material
interfaces and the solution depends on the choice for evaluating the
cell cross-sections and the leakage term.

![image](examples//02-general/02-reflected_circle/reflectedcircle.pdf)

#### reflected.mil {#reflected.mil .unnumbered .unnumbered}

This example shows the difference of the results obtained by the two
methods of evaluation cell cross sections, namely `XS_CENTER` and
`XS_MEAN`. By using commandline replacement, the same input solves a
reflected circle with both methods and then the thermal absorption cross
sections and the diffusion coefficients for each cell are shown. Even
though a nice "anti-aliasing" effect is obtained when using `XS_MEAN`,
the computational effort to compute the mean values may be many times
the effort needed to evaluate the parameters at the cell center without
giving a significant improvement in the solution. The usage of either
method has to be studied for each particular problem to be solved.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # reflected circle with different XS association
    # provide either XS_CENTER or XS_MEAN in the commandline

    PROBLEM DIMENSIONS 2 GROUPS 2 

    SCHEME $1

    x_bare_length = 100
    y_bare_length = x_bare_length
    # a low number of cells is given to study the cell XS
    # evaluation scheme
    x_cells = 40
    y_cells = x_cells

    # reflector width
    a = 20

    INCLUDE materials.mil

    ZONE fuel MATERIAL fuel     {
      X_CENTER x_bare_length/2
      Y_CENTER x_bare_length/2
      OUTER_RADIUS x_bare_length/2
    }

    ZONE refl MATERIAL reflector {
      X_CENTER x_bare_length/2
      Y_CENTER x_bare_length/2
      INNER_RADIUS x_bare_length/2-a
      OUTER_RADIUS x_bare_length/2
    }

    PRINT_FUNCTION SigmaT_2 D_1 flux_1 flux_2

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga reflected.mil XS_CENTER > xscenter.dat
    $ milonga reflected.mil XS_MEAN > xsmean.dat
    $ gnuplot xscenter.gnuplot
    $ gnuplot xsmean.gnuplot
    $ gnuplot 3dxscenter.gnuplot
    $ gnuplot 3dxsmean.gnuplot
    $

![image](examples//02-general/02-reflected_circle/sigmaacenter.pdf)

![image](examples//02-general/02-reflected_circle/Dcenter.pdf)

![image](examples//02-general/02-reflected_circle/centerfast.pdf)

![image](examples//02-general/02-reflected_circle/centerthermal.pdf)

![image](examples//02-general/02-reflected_circle/3dxscenter.pdf)

![image](examples//02-general/02-reflected_circle/sigmaamean.pdf)

![image](examples//02-general/02-reflected_circle/Dmean.pdf)

![image](examples//02-general/02-reflected_circle/meanfast.pdf)

![image](examples//02-general/02-reflected_circle/meanthermal.pdf)

![image](examples//02-general/02-reflected_circle/3dxsmean.pdf)

#### rods.mil {#rods.mil .unnumbered .unnumbered}

As with the slab, there are no symmetric reactors and thus full-core
calculations do not make sense. The following example includes some
absorbers that again do not coincide with cell borders, that would
represent most of the cases of interest. As with the slab, first the
continuous absorption cross section is given, showing the real position
of the absorbers. Afterward, the cell values using averages are given,
showing the antialiasing effect.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # reflected circle with some absorbers
    PROBLEM DIMENSIONS 2 GROUPS 2 

    x_bare_length = 100
    y_bare_length = x_bare_length
    x_cells = 40
    y_cells = x_cells

    a = 20 # reflector width

    INCLUDE materials.mil

    ZONE fuel MATERIAL fuel     {
      X_CENTER x_bare_length/2
      Y_CENTER x_bare_length/2
      OUTER_RADIUS x_bare_length/2
    }

    ZONE refl MATERIAL reflector {
      X_CENTER x_bare_length/2
      Y_CENTER x_bare_length/2
      INNER_RADIUS x_bare_length/2-a
      OUTER_RADIUS x_bare_length/2
    }

    # these zones do not coincide with cell borders on purpose
    ZONE rod1 MATERIAL rod       {
      X_MIN 51.5 X_MAX 72.2
      Y_MIN 20.1 Y_MAX 30.9
    }

    ZONE rod2 MATERIAL rod       {
      X_CENTER 65  Y_CENTER 70
      OUTER_RADIUS 14.75
    }

    ZONE rod3 MATERIAL rod       {
      X_CENTER 34.3
      Y_CENTER 57.4
      OUTER_RADIUS 4.12
    }

    # again, the cell XS do not need range
    PRINT_FUNCTION SigmaT_2 D_1 flux_1 flux_2

    # the continuous XS need ranges and steps. The first argument
    # is the energy group, the second is x and the third is y
    FILE cont cont.dat
    PRINT_FUNCTION FILE cont cSigmaT MIN 2 0 0 MAX 2 x_bare_length y_bare_length STEP 1 0.1*x_bare_length/x_cells 0.1*y_bare_length/x_cells

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga rods.mil > rods.dat
    $ gnuplot xsrods.gnuplot
    $ gnuplot 3drods.gnuplot
    $

![image](examples//02-general/02-reflected_circle/csigmaarods.pdf)

![image](examples//02-general/02-reflected_circle/sigmaarods.pdf)

![image](examples//02-general/02-reflected_circle/Drods.pdf)

![image](examples//02-general/02-reflected_circle/rodsfast.pdf)

![image](examples//02-general/02-reflected_circle/rodsthermal.pdf)

![image](examples//02-general/02-reflected_circle/3drods.pdf)

### IAEA 2D PWR Benchmark

This section shows how a classic IAEA benchmark for nuclear reactor
codes can be solved by milonga. Because the current version does not
work 100% with three-dimensional problems, the problem to be solved is
the so-called two-dimensional LWR Benchmark Problem with identification
11-A2 described in reference [@anl7416]. This problem represents the
mid-plane $z=190 \, \text{cm}$ of the 3D IAEA Benchmark Problem that
will be included as an example is future versions of milonga. The 3D
problem is a standard benchmark and its solution using with different
modern codes can be found, amongst other documents, in
[@imelda; @mosteller; @parcs3d]. The original reference [@anl7416]
provides also a collection of solutions obtained with legacy codes.

Quoting the original problem definition [@anl7416 page 437],

> Reduction of Source Situation:
>
> 1.  Two-group diffusion theory
>
> 2.  Two-dimensions $(x,y)$-geometry
>
> Two-Group Diffusion Equations:
>
> $$\begin{aligned}
> &-\nabla D_1 \nabla \phi_1 + (\Sigma_{a1} + \Sigma_{1 \rightarrow 2} + D_1 B_{z1}^2) \, \phi_1 = \frac{1}{\lambda} \, \nu \Sigma_{f2} \, \phi_2 \\
> &-\nabla D_2 \nabla \phi_2 + (\Sigma_{a2} + D_2 B_{z2}^2 ) \, \phi_2 = \Sigma_{1 \rightarrow 2} \phi_1\end{aligned}$$
>
> Data (see figure).
>
> Axial buckling $B_{zg}^2 = 0.8 \cdot 10^{-4}$ for all regions an
> energy groups
>
> Boundary conditions $J_g^{\text{in}} = 0$ no incoming current at
> external boundaries. For finite difference diffusion theory codes the
> following form is considered equivalent
>
> $$\frac{\partial \phi_g}{\partial n} = - \frac{0.4682}{D_g} \, \phi_g$$
> where $n$ is the outward directed normal to the surface. At symmetry
> boundaries:
>
> $$\frac{\partial \phi_g}{\partial n} = 0$$

![image](examples//02-general/03-iaea2dpwr/quarter.pdf)

![image](examples//02-general/03-iaea2dpwr/constants.pdf)

#### iaea2dpwr.mil {#iaea2dpwr.mil .unnumbered .unnumbered}

The following input file solves the 2D IAEA PWR benchmark described
above. The original problem asks for a number of results that can all be
computed by milonga, but for simplicity reasons only the multiplication
factor, flux 2D distribution $\phi(x,y)$ and radial profiles trough the
horizontal axis $\phi(x,0)$ and through the diagonal $\phi(x,x)$ as a
function of $x$ are given. The axial buckling is inserted as a
modification to the absorption cross section in the material file. The
diffusion equation is solved using the finite volumes scheme. As the
material interfaces coincide with cell borders, cell cross sections are
computed with the center method and currents are computed using the
diffusion coefficient evaluated at a distance $\epsilon$. An analysis of
how the results---including calculation times--change with the mesh size
can be easily performed by using parametric calculations over the
variable `c` that controls the mesh refinement, as shown in the next
section. Boundary conditions are set to null flux at the external
boundary and mirror conditions at $x=0$ and $y=0$ as to represent one
quarter of the core. The flux profile along the $x$ axis is computed
at $y=\Delta y/2$ instead of at $y=0$ to avoid spending time
interpolating by asking for the value the fluxes take exactly at a cell
center. The result is the same, because the mirror boundary condition
at $y=0$ sets  $\phi(x,-\Delta y/2)=\phi(x,\Delta y/2)$ and thus
$\phi(x,0)=\phi(x,\Delta y/2)$.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # IAEA 2D PWR benchmark
    # ANL-7416 Supplement 2, 1977, Argonne Code Center:
    # Benchmark Problem Book, page 437

    # two spatial dimensions and two energy groups
    PROBLEM DIMENSIONS 2 GROUPS 2

    # as the material interfaces coincide with cell borders,
    # the selection of XS_CENTER and XS_MEAN gives the same cell XS
    # also D_EPSILON and D_MEAN give the same currents
    # this options reduce computation time as they avoid integrations
    SCHEME XS_CENTER D_EPSILON

    # bare length in cm as defined in the problem
    x_bare_length = 170
    y_bare_length = 170

    # mesh refinement factor: the higher this value the finer the mesh
    c = 2.0

    x_cells = c*x_bare_length
    y_cells = x_cells
    deltax = x_bare_length/x_cells
    deltay = y_bare_length/y_cells

    # a quarter core in the first x-y quadrant
    BOUNDARY_CONDITIONS X_MIN MIRROR X_MAX NULL Y_MIN MIRROR Y_MAX NULL

    # zone defintion
    # new zones override previous ones when they overlap
    # note that Y_MIN defaults to zero
    ZONE refl_1    MATERIAL reflector  X_MIN   0 X_MAX  70 Y_MAX 170
    ZONE refl_2    MATERIAL reflector  X_MIN  70 X_MAX 110 Y_MAX 150
    ZONE refl_3    MATERIAL reflector  X_MIN 110 X_MAX 130 Y_MAX 130
    ZONE refl_4    MATERIAL reflector  X_MIN 130 X_MAX 150 Y_MAX 110
    ZONE refl_5    MATERIAL reflector  X_MIN 150 X_MAX 170 Y_MAX  70

    ZONE fuel2_1   MATERIAL fuel1      X_MIN   0 X_MAX  50 Y_MAX 150
    ZONE fuel2_2   MATERIAL fuel1      X_MIN  50 X_MAX  90 Y_MAX 130
    ZONE fuel2_3   MATERIAL fuel1      X_MIN  90 X_MAX 110 Y_MAX 110
    ZONE fuel2_4   MATERIAL fuel1      X_MIN 110 X_MAX 130 Y_MAX  90
    ZONE fuel2_5   MATERIAL fuel1      X_MIN 130 X_MAX 150 Y_MAX  50

    ZONE fuel1_1   MATERIAL fuel2      X_MIN   0 X_MAX  30 Y_MAX 130
    ZONE fuel1_2   MATERIAL fuel2      X_MIN  30 X_MAX  70 Y_MAX 110
    ZONE fuel1_3   MATERIAL fuel2      X_MIN  70 X_MAX 110 Y_MAX  70
    ZONE fuel1_4   MATERIAL fuel2      X_MIN 110 X_MAX 130 Y_MAX  30

    ZONE fuelrod_1 MATERIAL fuel2rod   X_MIN   0 Y_MIN   0 X_MAX  10 Y_MAX  10
    ZONE fuelrod_2 MATERIAL fuel2rod   X_MIN   0 Y_MIN  70 X_MAX  10 Y_MAX  90
    ZONE fuelrod_3 MATERIAL fuel2rod   X_MIN  70 Y_MIN   0 X_MAX  90 Y_MAX  10
    ZONE fuelrod_4 MATERIAL fuel2rod   X_MIN  70 Y_MIN  70 X_MAX  90 Y_MAX  90  

    # the materials are included in a different file
    # in order to be shared between other input files
    INCLUDE materials.mil

    # extract one-dimensional profiles from the fluxes
    # the deltay/2 is to save time in making
    # a multidimensional interpolation by asking
    # for values exactly at a cell center
    FUNCTION fastprofile(x) = flux_1(x, deltay/2)
    FUNCTION thermalprofile(x) = flux_2(x, deltay/2)

    FUNCTION fastdiagonal(x) = flux_1(x, x)
    FUNCTION thermaldiagonal(x) = flux_2(x, x)

    # output some information to the screen
    PRINT TEXT "keff = " keff FORMAT %.8lf
    PRINT TEXT "build_time = " build_time TEXT "seconds" FORMAT %.4lf
    PRINT TEXT "solve_time = " solve_time TEXT "seconds" FORMAT %.4lf
    PRINT TEXT "total_time = " build_time+solve_time TEXT "seconds" FORMAT %.4lf

    # two-dimensional flux distribution and total thermal XS
    FILE flux flux.dat
    PRINT_FUNCTION FILE flux flux_1 flux_2 SigmaT_2

    # flux profiles
    FILE profile profile.dat
    FILE diagonal diagonal.dat

    # these are algebraic functions, a range is mandatory
    # the range and step are selected again to avoid interpolation
    PRINT_FUNCTION FILE profile  fastprofile thermalprofile   MIN deltax/2 MAX x_bare_length-deltax/2 STEP deltax
    PRINT_FUNCTION FILE diagonal fastdiagonal thermaldiagonal MIN deltax/2 MAX x_bare_length-deltax/2 STEP deltax

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ cat materials.mil
    # materials for the IAEA PWR benchmark problem
    # a geometric bucking of 0.8e-4 in the z-direction
    # is taken into account by a D B_g^2 term in the
    # absorption cross section

    MATERIAL fuel1 {
      D_1         1.500
      D_2         0.400
      SigmaS_1.2  0.020
      SigmaA_1    0.010+1.5*0.8e-4
      SigmaA_2    0.080+0.4*0.8e-4
      nuSigmaF_2  0.135
    }

    MATERIAL fuel2 {
      D_1         1.500
      D_2         0.400
      SigmaS_1.2  0.020
      SigmaA_1    0.010+1.5*0.8e-4
      SigmaA_2    0.085+0.4*0.8e-4
      nuSigmaF_2  0.135
    }

    MATERIAL fuel2rod {
      D_1         1.500
      D_2         0.400
      SigmaS_1.2  0.020
      SigmaA_1    0.010+1.5*0.8e-4
      SigmaA_2    0.130+0.4*0.8e-4
      nuSigmaF_2  0.135
    }

    MATERIAL reflector {
      D_1         2.000
      D_2         0.300
      SigmaS_1.2  0.040+2.0*0.8e-4
      SigmaA_1    0.000+0.3*0.8e-4
      SigmaA_2    0.010
      nuSigmaF_2  0.000
    }
    $ milonga iaea2dpwr.mil
    keff = 	1.02990385	
    build_time = 	4.1996	seconds	
    solve_time = 	5.0606	seconds	
    total_time = 	9.2602	seconds	
    $ gnuplot iaea2d.gnuplot
    $ gnuplot profile.gnuplot
    $

\pagebreak
![image](examples//02-general/03-iaea2dpwr/sigmat2.pdf)

![image](examples//02-general/03-iaea2dpwr/fast.pdf)

![image](examples//02-general/03-iaea2dpwr/thermal.pdf)

![image](examples//02-general/03-iaea2dpwr/profile.pdf)

![image](examples//02-general/03-iaea2dpwr/diagonal.pdf)

#### circle.mil {#circle.mil .unnumbered .unnumbered}

A small variation of the 2D IAEA PWR Benchmark---that is also applicable
to the 3D problem---is the replacement of the external reflector
boundary by a circle instead of a coarse $17\times17$ rectangular
approximation. Of course, this requirement was not possible thirty five
years ago, but nowadays circular reflectors should be bread and butter
for core codes. This example shows one of the vectors of milonga's
design basis in action, namely the Independence of the problem geometry
and the spatial nodalization. A lower mesh refinement factor $c=0.5$ is
chosen to show how the circle is automatically translated into the
rectangular mesh. There is no need to use `XS_MEAN` because mainly the
only cells whose cross sections will differ from the center value are
those in the external boundary of the core and thus they have boundary
conditions equations and do not make use of the cell cross sections.
Note the reduction of the computation time because of the coarser mesh.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # Variation of the IAEA 2D PWR benchmark
    # quarter core with circular reflector

    PROBLEM DIMENSIONS 2 GROUPS 2
    SCHEME XS_CENTER D_EPSILON

    x_bare_length = 170
    y_bare_length = 170

    # mesh refinement factor
    c = 0.5

    x_cells = c*x_bare_length
    y_cells = x_cells
    deltax = x_bare_length/x_cells

    BOUNDARY_CONDITIONS X_MIN MIRROR X_MAX NULL Y_MIN MIRROR Y_MAX NULL

    # start with a circular reflector
    ZONE refl     MATERIAL reflector   X_CENTER 0 Y_CENTER 0 RADIUS 170

    # and then add the fuel elements in their correspnding square grids
    ZONE fuel2_1   MATERIAL fuel1      X_MIN   0 X_MAX  50 Y_MAX 150
    ZONE fuel2_2   MATERIAL fuel1      X_MIN  50 X_MAX  90 Y_MAX 130
    ZONE fuel2_3   MATERIAL fuel1      X_MIN  90 X_MAX 110 Y_MAX 110
    ZONE fuel2_4   MATERIAL fuel1      X_MIN 110 X_MAX 130 Y_MAX  90
    ZONE fuel2_5   MATERIAL fuel1      X_MIN 130 X_MAX 150 Y_MAX  50

    ZONE fuel1_1   MATERIAL fuel2      X_MIN   0 X_MAX  30 Y_MAX 130
    ZONE fuel1_2   MATERIAL fuel2      X_MIN  30 X_MAX  70 Y_MAX 110
    ZONE fuel1_3   MATERIAL fuel2      X_MIN  70 X_MAX 110 Y_MAX  70
    ZONE fuel1_4   MATERIAL fuel2      X_MIN 110 X_MAX 130 Y_MAX  30

    ZONE fuelrod_1 MATERIAL fuel2rod   X_MIN   0 Y_MIN   0 X_MAX  10 Y_MAX  10
    ZONE fuelrod_2 MATERIAL fuel2rod   X_MIN   0 Y_MIN  70 X_MAX  10 Y_MAX  90
    ZONE fuelrod_3 MATERIAL fuel2rod   X_MIN  70 Y_MIN   0 X_MAX  90 Y_MAX  10
    ZONE fuelrod_4 MATERIAL fuel2rod   X_MIN  70 Y_MIN  70 X_MAX  90 Y_MAX  90  

    INCLUDE materials.mil


    PRINT TEXT "keff = " keff FORMAT %.8lf
    PRINT TEXT "build_time = " build_time TEXT "seconds" FORMAT %.4lf
    PRINT TEXT "solve_time = " solve_time TEXT "seconds" FORMAT %.4lf
    PRINT TEXT "total_time = " build_time+solve_time TEXT "seconds" FORMAT %.4lf

    FILE flux circflux.dat
    PRINT_FUNCTION FILE flux flux_1 flux_2 SigmaT_2

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga circle.mil
    keff = 	1.03089876	
    build_time = 	0.2296	seconds	
    solve_time = 	0.1351	seconds	
    total_time = 	0.3647	seconds	
    $ gnuplot circle.gnuplot
    $

![image](examples//02-general/03-iaea2dpwr/circsigmat2.pdf)

![image](examples//02-general/03-iaea2dpwr/circfast.pdf)

![image](examples//02-general/03-iaea2dpwr/circthermal.pdf)

### Two-dimensional slice of a PHWR

The last example of the general-interest problems is a also
two-dimensional slice but this time of a fictitious Pressurized Heavy
Water Reactor. This solved case shows some other features of milonga
that may be handy when tackling similar problems. Although this example
is loosely based on a real reactor, the actual data is less than random
and of course the results presented here do not have any physical
significance.

In this kind of reactors, one of the most important parameters that
control the actual spatial cross-section dependance inside the core is
the fuel burnup, as the spatial distribution of this latter property
depends on the online fuel refueling and shuffling strategy. The core
consists of a number of vertical channels, each representing an
homogeneous mixture of fuel, coolant and moderator---that are separated
in PHWRs. The macroscopic cross sections depend on several parameters
such as temperatures, xenon and fuel burnup. In this problem, only
dependance on burnup is taken into account, whilst xenon and
temperatures dependance are illustrated in sections that are about to
come. There is a circular reflector surrounding the core and there are
four pipes for the moderator loop, whose temperature is not taken into
account. There are eighteen circular control rods, of which only six
circular are inserted at the considered plane. However, guide tubes for
the withdrawn rods are also considered. There are also fifteen lances
containing in-core instrumentation that introduce extra absorbing
materials in the core.

![image](examples//02-general/04-phwr/phwr.pdf)

#### 2dphwr.mil {#dphwr.mil .unnumbered .unnumbered}

First, the fuel burnup is read and interpolated from a text file that
gives the location of each channel and the burnup at the axial plane
considered. One-dimensional functions of burnup are read and
interpolated from a big multi-column file containing macroscopic cross
sections as a function of burnup. Materials are defined accordingly,
where burnup distribution is taken into account in the reactor core.
Then, the geometry is entered by using a mixture of absolute and
relative zones. Finally, some information is written to files
while $k_\text{eff}$ and solution time is written to the screen.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # 2D ficticious PHWR

    PROBLEM DIMENSIONS 2 GROUPS 2

    # because the burnup is not really interpolated but continuously
    # filled with the closest value (see below), it makes no sense
    # to waste CPU time using XS_MEAN 
    SCHEME XS_CENTER D_EPSILON

    # a bare square 7m long
    x_bare_length = 700
    y_bare_length = 700

    # square cells
    x_cells = 200
    y_cells = x_cells

    # the fuel burnup is given as a continuous function read from
    # a file, whose content can be seen in the terminal screen
    # interpolation method closest is selected to save computation time
    # and avoid extrapolation problems at the outer parts of the core
    # columns of burnup.dat are [cm] [cm] [MWday/ton]
    FUNCTION burnup(x,y) FILE burnup.dat INTERPOLATION closest

    # these functions give the core cross sections as a function of
    # the single parameter burnup (temperatures, xenon, etc are not
    # illustrated here but could be taken into account also)
    # XS units are 1/cm and bu units are MWday/ton
    FUNCTION D1(bu)         FILE xs.dat INTERPOLATION splines COLUMNS 1 2
    FUNCTION D2(bu)         FILE xs.dat INTERPOLATION splines COLUMNS 1 3
    FUNCTION SigmaA1(bu)    FILE xs.dat INTERPOLATION splines COLUMNS 1 4
    FUNCTION SigmaS12(bu)   FILE xs.dat INTERPOLATION splines COLUMNS 1 5
    FUNCTION SigmaS21(bu)   FILE xs.dat INTERPOLATION splines COLUMNS 1 6
    FUNCTION SigmaA2(bu)    FILE xs.dat INTERPOLATION splines COLUMNS 1 7
    FUNCTION nuSigmaF1(bu)  FILE xs.dat INTERPOLATION splines COLUMNS 1 8
    FUNCTION nuSigmaF2(bu)  FILE xs.dat INTERPOLATION splines COLUMNS 1 9

    # radial reflector properties are homogeneous
    MATERIAL reflector {
     D_1         1.405310E+01
     D_2         9.666628E-01
     SigmaA_1    4.107654E-06
     SigmaS_1.2  1.090063E-02
     SigmaS_2.1  1.088030E-04
     SigmaA_2    4.908041E-05
    }

    # core properties depend on burnup only
    MATERIAL fuelmod {
     D_1         D1(burnup(x,y))
     D_2         D2(burnup(x,y))
     SigmaA_1    SigmaA1(burnup(x,y))
     SigmaS_1.2  SigmaS12(burnup(x,y))
     SigmaS_2.1  SigmaS21(burnup(x,y))
     SigmaA_2    SigmaA2(burnup(x,y))
     nuSigmaF_1  nuSigmaF1(burnup(x,y))
     nuSigmaF_2  nuSigmaF2(burnup(x,y))
    }

    # incremental control rod XS
    MATERIAL rod {
      SigmaA_1      1.957196E-03
      SigmaA_2      7.339746E-03
      SigmaS_1.2   -7.986794E-04
      SigmaS_2.1    2.732369E-05
      nuSigmaF_1   -4.314411E-05
      nuSigmaF_2    5.649112E-04
    }

    # incremental control rod guide tube XS
    MATERIAL guide {
      SigmaA_1    -1.102107E-04
      SigmaA_2     9.838934E-05
      SigmaS_1.2  -1.378342E-04
      SigmaS_2.1   1.513809E-06
      nuSigmaF_1  -5.002725E-06
      nuSigmaF_2   1.110945E-05
    }

    # absoulte XS for the detectors
    MATERIAL detector {
      D_1           1
      D_2           1
      SigmaA_1      1e-3
      SigmaA_2      1e-2
      SigmaS_1.2    1e-4
      SigmaS_2.1    0
    }


    # first define the reflector (downcomer) as a big circle
    ZONE refl  MATERIAL reflector X_CENTER 350 Y_CENTER 350 RADIUS 345
    # filled with a smaller circle representing the core
    ZONE core  MATERIAL fuelmod   X_CENTER 350 Y_CENTER 350 RADIUS 305
    # four circles to hold the moderator pipes
    ZONE mod1  MATERIAL reflector X_CENTER 430 Y_CENTER 60  RADIUS 35
    ZONE mod2  MATERIAL reflector X_CENTER 640 Y_CENTER 270 RADIUS 35
    ZONE mod3  MATERIAL reflector X_CENTER 290 Y_CENTER 640 RADIUS 35
    ZONE mod4  MATERIAL reflector X_CENTER 60  Y_CENTER 450 RADIUS 35
     
    # guide tubes for the 18 control rods
    ZONE tube1  MATERIAL guide X_CENTER 217 Y_CENTER 389  RADIUS 15 INCREMENTAL
    ZONE tube2  MATERIAL guide X_CENTER 177 Y_CENTER 353  RADIUS 15 INCREMENTAL
    ZONE tube3  MATERIAL guide X_CENTER 272 Y_CENTER 377  RADIUS 15 INCREMENTAL
    ZONE tube4  MATERIAL guide X_CENTER 299 Y_CENTER 507  RADIUS 15 INCREMENTAL
    ZONE tube5  MATERIAL guide X_CENTER 245 Y_CENTER 471  RADIUS 15 INCREMENTAL
    ZONE tube6  MATERIAL guide X_CENTER 367 Y_CENTER 495  RADIUS 15 INCREMENTAL
    ZONE tube7  MATERIAL guide X_CENTER 435 Y_CENTER 436  RADIUS 15 INCREMENTAL
    ZONE tube8  MATERIAL guide X_CENTER 449 Y_CENTER 507  RADIUS 15 INCREMENTAL
    ZONE tube9  MATERIAL guide X_CENTER 408 Y_CENTER 412  RADIUS 15 INCREMENTAL
    ZONE tube10 MATERIAL guide X_CENTER 517 Y_CENTER 318  RADIUS 15 INCREMENTAL
    ZONE tube11 MATERIAL guide X_CENTER 503 Y_CENTER 389  RADIUS 15 INCREMENTAL
    ZONE tube12 MATERIAL guide X_CENTER 463 Y_CENTER 271  RADIUS 15 INCREMENTAL
    ZONE tube13 MATERIAL guide X_CENTER 381 Y_CENTER 224  RADIUS 15 INCREMENTAL
    ZONE tube14 MATERIAL guide X_CENTER 449 Y_CENTER 212  RADIUS 15 INCREMENTAL
    ZONE tube15 MATERIAL guide X_CENTER 381 Y_CENTER 271  RADIUS 15 INCREMENTAL
    ZONE tube16 MATERIAL guide X_CENTER 258 Y_CENTER 259  RADIUS 15 INCREMENTAL
    ZONE tube17 MATERIAL guide X_CENTER 258 Y_CENTER 200  RADIUS 15 INCREMENTAL
    ZONE tube18 MATERIAL guide X_CENTER 217 Y_CENTER 294  RADIUS 15 INCREMENTAL

    # six inserted control rods
    ZONE rod1 MATERIAL rod X_CENTER 177 Y_CENTER 353 RADIUS 15 INCREMENTAL
    ZONE rod2 MATERIAL rod X_CENTER 258 Y_CENTER 200 RADIUS 15 INCREMENTAL
    ZONE rod3 MATERIAL rod X_CENTER 449 Y_CENTER 212 RADIUS 15 INCREMENTAL
    ZONE rod4 MATERIAL rod X_CENTER 449 Y_CENTER 507 RADIUS 15 INCREMENTAL
    ZONE rod5 MATERIAL rod X_CENTER 245 Y_CENTER 471 RADIUS 15 INCREMENTAL
    ZONE rod6 MATERIAL rod X_CENTER 517 Y_CENTER 318 RADIUS 15 INCREMENTAL

    # fifteen in-core detectors
    ZONE det1  MATERIAL detector X_CENTER 299 Y_CENTER 365 RADIUS 5
    ZONE det2  MATERIAL detector X_CENTER 394 Y_CENTER 389 RADIUS 5
    ZONE det3  MATERIAL detector X_CENTER 367 Y_CENTER 342 RADIUS 5
    ZONE det4  MATERIAL detector X_CENTER 272 Y_CENTER 318 RADIUS 5
    ZONE det5  MATERIAL detector X_CENTER 367 Y_CENTER 459 RADIUS 5
    ZONE det6  MATERIAL detector X_CENTER 435 Y_CENTER 294 RADIUS 5
    ZONE det7  MATERIAL detector X_CENTER 231 Y_CENTER 412 RADIUS 5
    ZONE det8  MATERIAL detector X_CENTER 463 Y_CENTER 436 RADIUS 5
    ZONE det9  MATERIAL detector X_CENTER 326 Y_CENTER 224 RADIUS 5
    ZONE det10 MATERIAL detector X_CENTER 177 Y_CENTER 483 RADIUS 5
    ZONE det11 MATERIAL detector X_CENTER 558 Y_CENTER 459 RADIUS 5
    ZONE det12 MATERIAL detector X_CENTER 340 Y_CENTER 129 RADIUS 5
    ZONE det13 MATERIAL detector X_CENTER 190 Y_CENTER 247 RADIUS 5
    ZONE det14 MATERIAL detector X_CENTER 367 Y_CENTER 554 RADIUS 5
    ZONE det15 MATERIAL detector X_CENTER 531 Y_CENTER 271 RADIUS 5


    # output files
    FILE xsburn xsburnup.dat            # 1D-dependance of XS with burnup
    FILE burnup burnup-interpolated.dat # 2D distribution of interpolated burnup
    FILE xsdist xsdist.dat              # 2D distribution of cells XS
    FILE solution solution.dat          # 2D flux distribution

    # althought range is optional, it is given to explicitly show that
    # these functions are interpolated as a function of the burnup
    # that may take virtually any range between 0 and 10000
    # in the gnuplot input, file xs.dat with the actual definition
    # points is plotted over the continuous spline interpolation
    PRINT_FUNCTION FILE xsburn SigmaA1 SigmaA2 nuSigmaF1 nuSigmaF2 SigmaS12 SigmaS21 MIN 0 MAX 10000 STEP 50

    # burnup distribution as interpolated from burnup.dat
    # range is optional but given in order to explicitly show how
    # the interpolation method "closest" work
    PRINT_FUNCTION FILE burnup burnup MIN 0 0 MAX x_bare_length y_bare_length STEP x_bare_length/x_cells  y_bare_length/y_cells

    # cell XS distribution
    PRINT_FUNCTION FILE xsdist nuSigmaF_1 nuSigmaF_2 SigmaA_1 SigmaA_2

    # flux distribution
    PRINT_FUNCTION FILE solution flux_1 flux_2 


    PRINT TEXT "keff = " keff FORMAT %.6lf
    PRINT TEXT "time = " build_time TEXT " + " solve_time TEXT " = " build_time+solve_time TEXT "seconds" FORMAT %.2lf SEPARATOR " "

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ head burnup.dat
    # the first columns have the location of the channels centers in
    # the x-y plane, and the third one has the burnup in MW-day/ton
    354	353	5144.0
    381	353	9069.8
    367	377	6109.6
    340	377	7906.5
    326	353	6500.9
    340	330	8620.4
    367	330	7354.4
    394	377	4259.6
    $ milonga 2dphwr.mil
    keff = 	1.003979	
    time =  36.07  +  1.77  =  37.84 seconds 
    $ gnuplot xs.gnuplot
    $ gnuplot 2dxs.gnuplot
    $ gnuplot solution.gnuplot
    $

![image](examples//02-general/04-phwr/absxsbu.pdf)

![image](examples//02-general/04-phwr/fisxsbu.pdf)

![image](examples//02-general/04-phwr/scatxsbu.pdf)

![image](examples//02-general/04-phwr/burnup.pdf)

![image](examples//02-general/04-phwr/fisdist.pdf)

![image](examples//02-general/04-phwr/absdist.pdf)

![image](examples//02-general/04-phwr/fast.pdf)

![image](examples//02-general/04-phwr/thermal.pdf)

Parametric problems
-------------------

This section focuses on parametric calculations, i.e. systematically
varying one or more parameters to obtain results as a function of them.
This feature may be used for example to perform sensitivity analysis,
build design maps and optimize parameters. The problems shown below have
analytical solutions and the parameter varied systematically has to do
the geometry and the nodalization. This way, the problems are kept
simple, the expected solution is known *a priori* and some insights
about how milonga works and what results are to be expected are shown.
Engineering analysis would include higher dimensions and a few groups of
energy making use of parametric studies on compositions, temperatures,
burnable poisons, etc. Actually, at least one PhD Thesis should be
written using this tool. The more, the merrier.

### Grid size effects

The three examples that follow show solve again the one-speed bare slab,
bare square and a bare three-dimensional cube, but this time the
numerical multiplication factor is compared to the
analytical $k_\text{eff}$ as a function of the number of cells used in
the spatial discretization. All the cases use the same material---an
hypothetical fuel mixture with $k_\infty=1.2$---spanning the whole bare
space. In the last example, a comparison between the results obtained
for the three dimensions as a function of the problem size is performed.

#### 1d.mil {#d.mil .unnumbered .unnumbered}

A one-speed one-dimensional bare homogeneous slab is solved by
parametrically varying the number of cells. The file `materials.mil` is
shared between the three examples of this section, so it is shown in the
terminal view. For each step, the output is a single line in the
standard output. It consists of eight columns, namely the number of
cells in the $x$-direction, the size of the problem (equal to the number
of cells in this case but different for higher-order problems), the
numerical $k_\text{eff}$, the analytical $k_\text{eff}$ (that does not
depend on $\texttt{x\_cells}$ but is printed in each step to ease the
plot procedure), the time needed to build the matrices, the time needed
to solve the eigenvalue problem and the total time. Because no
`MAX_DAUGHTERS` keyword was entered, the calculation is done in series,
i.e. each step starts only after the last one has finished. Therefore,
the output lines will be ordered by increasing `x_cells`. If this had
not be the case, output should have been written to files and then
ordered from the shell, as explained in
chapter [3](#cap:input){reference-type="ref" reference="cap:input"}.
Conclusions about the results shown in the figures are up to the user,
but note that even though the problem size scales linearly, square
matrices of size $N\times N$ usually scale as $N^2$.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study of a bare slab with respect to
    # the number of cells
    PROBLEM DIMENSIONS 1 GROUPS 1

    # perform a parametric study on x_cells from 10
    # up to 5000 with an additive increment of 50
    PARAMETRIC x_cells 10 5000 50

    SCHEME DIFFERENCES XS_CENTER GRAD_D_LOCAL

    # 100 arbitrary units width
    x_bare_length = 100
    # this assignment is ignored in a parametric
    # calculation, however if one wants to comment
    # the PARAMETRIC keyword above, it is handy to
    # maintain this line to prevent milonga from failing
    x_cells = 100

    # materials and zones
    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel

    # analytical solution (variables are defined in
    # the materials.mil file)
    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + D*(pi/x_bare_length)^2)

    # in each step print one line as the output to
    # obtain how all these stuff depend on x_cells
    PRINT x_cells x_cells keff k keff-k build_time solve_time build_time+solve_time

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ cat materials.mil
    D = 1
    SigmaA = 1e-2
    nuSigmaF = 1.2e-2

    MATERIAL fuel {
      D_1        D
      SigmaA_1   SigmaA
      nuSigmaF_1 nuSigmaF
    }
    $ milonga 1d.mil > 1d.dat
    $ gnuplot 1d.gnuplot
    $

![image](examples//03-parametric/01-grid/1dk.pdf)

![image](examples//03-parametric/01-grid/1derror.pdf)

![image](examples//03-parametric/01-grid/1dtime.pdf)

#### 2d.mil {#d.mil-1 .unnumbered .unnumbered}

The two-dimensional bare slab is studied parametrically as a function of
the number of cells in the $x$-direction. In this case, the number of
cells in the $y$-direction is the same as in the $x$-direction, so the
problem size scales as the square of `x_cells`. Same figure as for the
one-dimensional cases are produced. For fun, compare how the results
depend on the number of cells with respect to the previous case, taking
into account the same remark regarding $N\times N$ and $N^2$.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study of a bare square with respect to
    # the number of cells
    PROBLEM DIMENSIONS 2 GROUPS 1

    # the problem size scales with x_cells^2 so
    # the range has to be chosen carefully to avoid
    # running out of memory
    PARAMETRIC x_cells 10 250 10

    SCHEME DIFFERENCES XS_CENTER GRAD_D_LOCAL

    x_bare_length = 100
    y_bare_length = x_bare_length

    # these two lines are ignored in parametric mode,
    # but again it is handy to keep them
    x_cells = 100
    y_cells = x_cells

    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel

    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + 2*D*(pi/x_bare_length)^2)

    PRINT x_cells x_cells^2 keff k keff-k build_time solve_time build_time+solve_time

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga 2d.mil > 2d.dat
    $ gnuplot 2d.gnuplot
    $

![image](examples//03-parametric/01-grid/2dk.pdf)

![image](examples//03-parametric/01-grid/2derror.pdf)

![image](examples//03-parametric/01-grid/2dtime.pdf)

#### 3d.mil {#d.mil-2 .unnumbered .unnumbered}

Lastly, the first three-dimensional example appears. It is a bare
one-dimensional box, whose effective multiplication factor is computed
for different spatial nodalizations. In each direction the number of
cells are equal to `x_cells`. After the usual three figures with the
$k_\text{eff}$, the error and the solution times, two figures summing up
the situation for each example of the section are provided. It can be
seen how, for the same matrix size lower dimensional problems are faster
and more accurate, as expected.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study of a bare cube with respect to
    # the number of cells
    PROBLEM DIMENSIONS 3 GROUPS 1

    PARAMETRIC x_cells 5 40 1

    SCHEME DIFFERENCES XS_CENTER GRAD_D_LOCAL

    x_bare_length = 100
    y_bare_length = x_bare_length
    z_bare_length = x_bare_length

    x_cells = 25
    y_cells = x_cells
    z_cells = x_cells

    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel

    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + 3*D*(pi/x_bare_length)^2)

    PRINT x_cells x_cells^3 keff k keff-k build_time solve_time build_time+solve_time

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga 3d.mil > 3d.dat
    $ gnuplot 3d.gnuplot
    $ gnuplot dims.gnuplot
    $

![image](examples//03-parametric/01-grid/3dk.pdf)

![image](examples//03-parametric/01-grid/3derror.pdf)

![image](examples//03-parametric/01-grid/3dtime.pdf)

![image](examples//03-parametric/01-grid/errors.pdf)

![image](examples//03-parametric/01-grid/times.pdf)

#### 3dinverse.mil {#dinverse.mil .unnumbered .unnumbered}

As a bonus track, another three-dimensional example is included. This
time, an illustration of how knowledge about numerical methods,
preconditioners and over-relaxation parameters can drastically reduce
computation times. This problem passes run-time options to SLEPc in the
commandline to select the Krylov-subspace solver GMRES, a SOR
preconditioner and to set some other optional parameters, resulting in a
reduction of the problem time for big problems when compared to the
previous case, as can be seen in the output figure. The SLEPc parameters
shown in the terminal view are not important *per se*. This example is
about the influence they have in milonga's performance. The actual
meaning of the options and its impact should be addressed on the one
hand by reading SLEPc's documentation [@slepc-users-manual] and on the
other hand understanding how milonga builds the associated matrices, as
explained in chapter [2](#cap:equations){reference-type="ref"
reference="cap:equations"} of this document.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study of a bare cube with respect to
    # the number of cells using the "inverse k" formulation
    # with sor preconditioner and gmres solver
    PROBLEM DIMENSIONS 3 GROUPS 1

    PARAMETRIC x_cells 5 40 1

    # for the SLEPc options in commandline to work, we have to solve
    # the problem were the eigenvalue is equal to 1/keff
    SOLVER EIGENVALUE_INVERSE_K
    SCHEME DIFFERENCES XS_CENTER GRAD_D_LOCAL

    x_bare_length = 100
    y_bare_length = x_bare_length
    z_bare_length = x_bare_length

    x_cells = 25
    y_cells = x_cells
    z_cells = x_cells

    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel

    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + 3*D*(pi/x_bare_length)^2)

    PRINT x_cells x_cells^3 keff k keff-k build_time solve_time build_time+solve_time

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga 3dinverse.mil -st_ksp_type gmres -st_pc_type sor -st_type sinvert -st_pc_sor_omega 1.7 -eps_ncv 7 > 3dinverse.dat
    $ gnuplot 3dinverse.gnuplot
    $

![image](examples//03-parametric/01-grid/3dinverse.pdf)

### Discrete boundary conditions effects

The following example presents another example of parametric
calculation. Again a bare slab is solved but this time, the number of
cells and the bare length are kept constant while the width $a$ is
changed parametrically. This leads to the general situation where the
external boundary does not coincide with cell borders, as shown in the
figure. The numerical multiplicative factor is again compared to the
theoretical one. The latter depends continuously on the slab width while
the second does not because of the discrete boundary condition
equations. The output figure shows how the error changes with the slab
width. Keep in mind that the cells are two units width and refer to
figures [\[fig:bc-planar\]](#fig:bc-planar){reference-type="ref"
reference="fig:bc-planar"}
and [\[fig:bc-planarneg\]](#fig:bc-planarneg){reference-type="ref"
reference="fig:bc-planarneg"} when analyzing the result.

![image](examples//03-parametric/02-bc/bc-fig.pdf)

#### bc.mil {#bc.mil .unnumbered .unnumbered}

This is a rather simple input. Zone `fuel` is defined to be in the range
$[0,a]$, while `x_bare_length` and `x_cells` are kept constant.
Variable `a` is changed parametrically and the error as a function
of $a$ is written to the standard output.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # effect of the application of discrete boundary
    # conditions on planar surfaces
    PROBLEM DIMENSIONS 1 GROUPS 1

    # sweep the slab width from 80 to 100
    PARAMETRIC a 80 100 0.1

    x_bare_length = 100
    # 50 cells give cells of width 2
    x_cells = 50

    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel X_MIN 0 X_MAX a

    pi = 4*atan(1)
    k = nuSigmaF/(SigmaA + D*(pi/a)^2)

    PRINT a keff k keff-k

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga bc.mil > bc.dat
    $ gnuplot bc.gnuplot
    $

![image](examples//03-parametric/02-bc/bc.pdf)

### Circle quadrature

In a similar fashion, this section shows what happens when a continuous
geometry does not coincide with a discrete mesh. This time, a bare
circle with a varying radius $r$ is solved and the multiplicative factor
vs. the radius is shown. Because of the complexity the many possible
combinations of interferences between a circular domain and a square
mesh, the dependance of $k_\text{eff}$ with $r$ is non-trivial. The
first example uses a $50\times 50$ mesh and the second a
$200 \times 200$, to show that even though there are ripples in the
numerical solution, it should converge to the analytical one for
infinite number of cells.

![image](examples//03-parametric/03-circle/quadrature.pdf)

#### circle50.mil {#circle50.mil .unnumbered .unnumbered}

The input below solved a bare circle of radius $r$ immersed in a square
bare domain of size $100 \times 100$ discretized with a $50 \times 50$
mesh. The values of the numerical an analytical $k_\text{eff}$ along
with the relative error are plotted as a function of the radius $r$.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study on keff vs radius
    # in a bare circle with a square mesh
    PROBLEM DIMENSIONS 2 GROUPS 1

    PARAMETRIC r 30 40 0.02

    x_bare_length = 100
    y_bare_length = x_bare_length

    x_cells = 50
    y_cells = x_cells

    INCLUDE materials.mil

    # centered circle of radius r
    ZONE fuel MATERIAL fuel X_CENTER x_bare_length/2 Y_CENTER y_bare_length/2 RADIUS r

    k = nuSigmaF/(SigmaA + D*(2.40483/r)^2)

    PRINT r keff k (keff-k)/k

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga circle50.mil > circle50.dat
    $ gnuplot circle50.gnuplot
    $

![image](examples//03-parametric/03-circle/keff50.pdf)

![image](examples//03-parametric/03-circle/error50.pdf)

#### circle200.mil {#circle200.mil .unnumbered .unnumbered}

Essentially the same problem is solved using a $200 \times 200$ mesh. A
comparison between the solution found with this input with respect to
the previous example is provided. Only a small interval is studied
because of running time reasons.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # parametric study on keff vs radius
    # in a bare circle with a square mesh
    PROBLEM DIMENSIONS 2 GROUPS 1

    PARAMETRIC r 34 36 0.02

    x_bare_length = 100
    y_bare_length = x_bare_length

    x_cells = 200
    y_cells = x_cells

    INCLUDE materials.mil
    ZONE fuel MATERIAL fuel X_CENTER x_bare_length/2 Y_CENTER y_bare_length/2 RADIUS r

    k = nuSigmaF/(SigmaA + D*(2.40483/r)^2)

    PRINT r keff k (keff-k)/k

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga circle200.mil > circle200.dat
    $ gnuplot circle200.gnuplot
    $

![image](examples//03-parametric/03-circle/keff200.pdf)

![image](examples//03-parametric/03-circle/error200.pdf)

### Control rod cusp problem

This example is an extension to the two-zone slab introduced in
section [4.1.6](#ex:twozone){reference-type="ref"
reference="ex:twozone"}. Now, the width $a$ of the absorber material
with $k_\infty<1$ is varied from 20 to 80 with and increment of
10$^{-2}$. The effective multiplication factor as computed from the
numerical problem is compared to the analytical $k_\text{eff}$ given by
the critical condition

$$\begin{aligned}
 &\sqrt{ D_\text{I} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right)} \cdot \tan \left[ \sqrt{ \frac{1}{D_\text{II}} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}}\right) }  \cdot \left( b - a \right) \right] + \\
 &\quad \quad \quad \sqrt{ D_\text{II} \left( \Sigma_{a\text{II}} - \frac{\nu\Sigma_{f\text{II}}}{k_\text{eff}} \right) } \cdot \tanh \left[ \sqrt{ \frac{1}{D_\text{I}} \left( \frac{\nu\Sigma_{f\text{I}}}{k_\text{eff}} - \Sigma_{a\text{I}} \right) } \cdot a \right] = 0\end{aligned}$$

Of course, $k_\text{eff}(a)$ should be a continuous and smooth function
of $a$. Indeed, the analytical solution is. However, depending on the
number of cells, the numerical scheme and the method of cell
cross-section association, the numerical solution may present cusps or
even discontinuities with respect to $a$. The user is encouraged to try
different nodalizations and options to analyze the results, sharing her
findings with milonga's author. Actually, the problem of cross-section
dilution was already studied using SLEPc [@Gilbert].

#### cusp.mil {#cusp.mil .unnumbered .unnumbered}

The number of cells selected is rather small to stress the effects of
the cusps on $k_\text{eff}(a)$ because of the dilution of cross sections
into finite cells.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # study of the effects of the rod cusp problem
    # using milonga's parametric capabilities

    PROBLEM DIMENSIONS 1 GROUPS 1

    # parameter a is varied from 20 to 80 with a 0.01 step
    PARAMETRIC a 20 80 0.01

    # this assignment is ignore when making a parametric
    # calculation on a
    a = 50
    b = 100
    x_bare_length = b

    # a small number of cells is selected in order
    # to stress the cusp problem
    x_cells = 25

    # scheme options
    ##################################
    # the user is urged to test all these combinations
    # and analyze the results obtained according to
    # the explanation given in milonga's documentation
    SCHEME VOLUMES XS_MEAN D_MEAN
    # SCHEME VOLUMES XS_MEAN D_EPSILON
    # SCHEME VOLUMES XS_CENTER D_MEAN
    # SCHEME VOLUMES XS_CENTER D_EPSILON
    # SCHEME DIFFERENCES XS_MEAN GRAD_D_NEIGHBORS
    # SCHEME DIFFERENCES XS_MEAN GRAD_D_LOCAL
    # SCHEME DIFFERENCES XS_CENTER GRAD_D_NEIGHBORS
    # SCHEME DIFFERENCES XS_CENTER GRAD_D_LOCAL


    # cross sections
    ##################################
    nu1 = 0.010
    nu2 = 0.015
    a1 = 0.015
    a2 = 0.010
    d1 = 0.5   
    d2 = 1.5

    MATERIAL fuel {
      D_1         d2
      SigmaT_1    a2
      nuSigmaF_1  nu2
    }

    MATERIAL rod {
      D_1         d1
      SigmaT_1    a1
      nuSigmaF_1  nu1
    }

    # zones
    ##################################
    ZONE rod  MATERIAL rod    X_MIN 0   X_MAX a
    ZONE fuel MATERIAL fuel   X_MIN a   X_MAX b

    # analytical multiplication factor
    ##################################
    # interval to look for the analytical k
    eps = 0.02 + 0.02*(a-20)/100

    # evaluate k from the analytical critical condition
    # taking as an initial guess the numerical keff
    k = keff
    k = root(sqrt(d1*abs(nu1/k-a1)) * tan(sqrt((1/d2)*abs(a2-nu2/k))*(b-a)) + sqrt(d2*abs(a2-nu2/k)) * tanh(sqrt((1/d1)*abs(nu1/k-a1))*a), k, keff-eps, keff+eps)

    PRINT a keff k

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga cusp.mil > cusp.dat
    $ gnuplot cusp.gnuplot
    $

![image](examples//03-parametric/04-cusp/keff.pdf)

![image](examples//03-parametric/04-cusp/zoom.pdf)

![image](examples//03-parametric/04-cusp/error.pdf)

Non-linear problems
-------------------

As stated in page , non-linear problems arise when the coefficients of
the diffusion equation depend themselves on the flux. To solve this kind
of problems, an iterative scheme may be used in the hope that it will
converge to the desired solution after a number of steps. The two
examples given in this sections deal with non-linear problems. The first
one illustrates how milonga treats xenon as a neutronic poison and the
second one implements a simple algorithm to find the position of a
control rod that renders a reactor critical. Again, these examples are
just an illustration of milonga features and do not posses any
particular engineering content whatsoever.

### Xenon effects

This example shows how to use milonga's capabilities regarding xenon
poisoning by solving again our good old friend, the bare one-speed
one-meter width slab. An iterative scheme is used, where the actual flux
is computed using the steady-state $^{135}$Xe distribution computed in
the previous step. A a power setpoint and a non-zero value
of $E \Sigma_f$ for at least one material has to be entered in order for
milonga to be able to compute a dimensional value for the flux
distribution $\phi(\ensuremath\mathbf{r})$ and therefore evaluate the
xenon distribution this flux distribution gives. Also, function
`xenon(x)` should be used when giving cross sections dependance to
define the xenon feedback effects. Calculation may end after a fixed
number of steps---as the case illustrated---or by setting a flag when
two successive $k_\text{eff}$'s differ less than a certain threshold.

#### xenon.mil {#xenon.mil .unnumbered .unnumbered}

A bare slab of one meter width is solved taking into account xenon
feedback effects. Units of `xenon(x)` are inverse cubic length, in
whatever units `x_bare_length` is. That is to say, this function gives a
density of $^{135}$Xe nuclei per volume. Multiplicative coefficients
should have appropriate units in order to recover inverse length for the
macroscopic cross sections. On the other hand, $E\Sigma_f$ has units of
energy times inverse length, where energy units should be compatible
with the units of the power setpoint. Default values for $\nu$, yields
and microscopic $^{135}$Xe neutron absorption are used. The
multiplicative factor vs. the step number---the first step considers no
poisoning---and the flux and xenon distribution for the first seven
steps are shown as outputs. Note that the xenon distribution for each
step is the resulting distribution and not the one used when computing
the flux. Again, for step one the xenon distribution used for
computing $\phi(x)$ is assumed to be identically zero.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # example illustrating how to include xenon
    # feedback effects in milonga

    # simple one-speed slab problem with a
    # small number of cells and center XS
    PROBLEM DIMENSIONS 1 GROUPS 1
    SCHEME VOLUMES XS_CENTER D_EPSILON
    x_bare_length = 100   # in cm
    x_cells = 25

    # file out is defined as step, so actually
    # one file per step will be created, called
    # out.n where n is the step number
    FILE out out STEP

    # we will perform fifteen steps
    # alternatively, a convergence check can be made by
    # comparing two sucessive keff and setting the variable
    # done to true
    static_iterations = 15


    # power setpoint: being a one dimensional problem, it
    # is a density per unit area perpendicular to x axis
    # length are in cm, so for a cubic reactor of size
    # 100x100x100 with a total power of 1 MW with energy
    # in joules and time in seconds, the power density
    # setpoint should be 
    # 1e6/(100x100) = 100
    # the power setpoint is mandatory for xenon calculations
    power = 100

    # keep the computed keff in the last step (before the
    # first ZONE keyword) to implement a convergence check
    kold = keff

    # a single zone spanning the whole bare length
    ZONE fuel MATERIAL fuel 

    # write into the output file flux, xenon concentration an XS
    PRINT_FUNCTION FILE out flux_1 xenon nuSigmaF_1 SigmaA_1
    # write into standard output the step, actual keff and
    # difference with respect to the last one
    PRINT static_step keff keff-kold

    # uncoment to quit whenever convergence is attained
    # instead of when reaching static_iterations
    # done = less(abs(keff-kold), 1e-5)

    # function xenon(x) contains atoms per volumetric units
    # ESigmaF has units of joules/cm
    MATERIAL fuel {
      D_1         "1       + 5e-18*xenon(x)"
      SigmaA_1    "4e-3    + 3e-17*xenon(x)"
      nuSigmaF_1  "nu*4e-3 - 5e-17*xenon(x)"
      ESigmaF_1   "200e6*1.6e-19*(2e-3 - 1e-17*xenon(x)/nu)"
    }

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga xenon.mil > xe.dat
    $ gnuplot xenon.gnuplot
    $

![image](examples//04-nonlinear/01-xenon/steps.pdf)

![image](examples//04-nonlinear/01-xenon/fluxes.pdf)

![image](examples//04-nonlinear/01-xenon/xenon.pdf)

### Criticallity with a control rod

The next example shows another kind of non-linear diffusion problem. It
consists of locating the position of a control rod such that the
resulting ensemble is as critical as possible. Again, this is just an
illustration of milonga's capabilities and not a real application. The
problem is a two-group two-dimensional horizontal core with a radial
reflector and four control rods entering from above. Two control rods
are fixed in space and the insertion of the other two is controlled by a
variable that is modified in each iteration according to the sign of the
resulting static reactivity.

#### rod.mil {#rod.mil-1 .unnumbered .unnumbered}

An initial position of 50% insertion is given by setting the initial
condition for variable `pos` to 0.5 by using the postfix `_init`. Then,
after solving the eigenvalue problem with the zone corresponding to the
two control rods defined in terms of the `pos` variable, the insertion
is updated as

$$\text{\texttt{pos}} \leftarrow \text{\texttt{pos}} + 10 \cdot (k_\text{eff}-1)$$
and a new iteration is performed. The flux distribution is written to
file `flux.dat`, that is not defined as `STEP` but the print instruction
has a condition that writes the information only if `static_step` is
equal to `static_iterations`, i.e. only in the last step.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # example of a simple algorithm to make a reactor
    # critical by moving a control rod

    PROBLEM DIMENSIONS 2 GROUPS 2

    static_iterations = 10 # number of static steps

    # bare lengths and nodalization
    x_bare_length = 600
    y_bare_length = 650
    x_cells = 100
    y_cells = 108

    # initial value (only for static_step = 1) of variable pos,
    # that is implicitly defined by this assignment also
    pos_init = 0.5

    MATERIAL fuel {
      D_1         1.500  D_2         0.400
      SigmaS_1.2  0.020  SigmaA_1    0.010
      SigmaA_2    0.085  nuSigmaF_2  0.130
    }

    MATERIAL fuelrod {
      D_1         1.500  D_2         0.400
      SigmaS_1.2  0.020  SigmaA_1    0.010
      SigmaA_2    0.130  nuSigmaF_2  0.130
    }

    MATERIAL reflector {
      D_1         2.000  D_2         0.300
      SigmaS_1.2  0.040  SigmaA_1    0.000
      SigmaA_2    0.010  nuSigmaF_2  0.000
    }

    # geometry definition
    # outer reflector
    ZONE refl MATERIAL reflector X_CENTER 300 Y_CENTER 300  OUTER_RADIUS 300 INNER_RADIUS 250
    # reactor core
    ZONE fuel MATERIAL fuel      X_CENTER 300 Y_CENTER 300  OUTER_RADIUS 250

    # two control rods whose insertion is defined by variable pos
    ZONE rod1 MATERIAL fuelrod   X_MIN 200 X_MAX 240 Y_MAX 700 Y_MIN 600*(1-pos)
    ZONE rod2 MATERIAL fuelrod   X_MIN 300 X_MAX 340 Y_MAX 700 Y_MIN 600*(1-pos)

    # two fixed control rods
    ZONE rod3 MATERIAL fuelrod   X_MIN 400 X_MAX 440 Y_MAX 700 Y_MIN 250
    ZONE rod4 MATERIAL fuelrod   X_MIN 150 X_MAX 160 Y_MAX 700 Y_MIN 150

    # write to standard output how the reactivity evolves as
    # the control rods are moved step by step
    PRINT static_step pos keff  (keff-1)/keff*1e5

    # file that will contain the last (critical) flux distribution
    # note that it is not defined with the STEP keyword as the PRINT_FUNCTION
    # instruction below makes sure only the last flux distribution is printed
    FILE flux flux.dat

    # the PRINT_FUNCTION is executed only if the condition evaluates
    # to true, i.e  only at  the last step of the iteration
    PRINT_FUNCTION CONDITION equal(static_step,static_iterations) FILE flux flux_1 flux_2 SigmaA_2 

    # algorithm to update the position of the control rod according
    # to the sign of the last computed reactivity
    pos = pos + 10*(keff-1)

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ milonga rod.mil
    1.000000e+00	5.000000e-01	1.007512e+00	7.456008e+02	
    2.000000e+00	5.000000e-01	1.007512e+00	7.456008e+02	
    3.000000e+00	5.751202e-01	1.005309e+00	5.281038e+02	
    4.000000e+00	6.282109e-01	1.002804e+00	2.796322e+02	
    5.000000e+00	6.562526e-01	1.000962e+00	9.608342e+01	
    6.000000e+00	6.658702e-01	1.000259e+00	2.592735e+01	
    7.000000e+00	6.684636e-01	1.000134e+00	1.337936e+01	
    8.000000e+00	6.698017e-01	1.000013e+00	1.254712e+00	
    9.000000e+00	6.699271e-01	9.999950e-01	-4.995780e-01	
    1.000000e+01	6.698772e-01	9.999950e-01	-4.995780e-01	
    $ gnuplot rod.gnuplot
    $

![image](examples//04-nonlinear/02-rod/sigmaa.pdf)

![image](examples//04-nonlinear/02-rod/flux.pdf)

Coupled calculations
--------------------

A coupled calculation is needed whenever there is at least one parameter
that influences the solution of the diffusion equation but is outside
milonga's scope. Therefore, an iterative scheme should be designed where
information between milonga and one or more external codes is to be
exchanged. For sure, coupled calculations are one of the most complex
task a nuclear engineer may have to perform and, as such, this is the
last section of the examples chapter. In these cases, the input files
may get rather cumbersome.


There are a number of ways of passing information from one process to
another, being the most straightforward---but also the most
inefficient---writing and reading from plain files. Modern operating
systems provide inter-process communication capabilities through kernel
system calls that can be used to couple calculation codes running in the
same host. Distributed coupled calculations may use network connections
to exchange information.

Although some file-based coupling may be implemented, milonga's coupling
scheme---actually provided by the wasora framework---is mainly based on
POSIX shared memory objects and semaphores. The main objective of this
scheme is to couple codes belonging to the wasora suite. Nevertheless,
coupled calculations with ad-hoc codes is also possible. Remote coupling
based on TCP connections will be implemented in future versions.

Anyhow, presenting examples about coupled calculations is a rather
difficult task as not only should milonga's features and characteristics
be discussed, but also the external program computing and coupling
capabilities should be explained. This section shows one single example
regarding a coupled calculation with thermalhydraulic feedback using
RELAP [@RELAP5VOL1] as the plant code.

### 1D core coupled with RELAP

This example solves a very simplified problem that involves the
computation of power, fuel and coolant temperature distributions in a
simplified model of a nuclear reactor core. A two-group one-dimensional
neutronic model computed by milonga is coupled to a one-dimensional
thermalhydraulic RELAP model computed. RELAP is a standard code
developed by the US N.R.C. to evaluate power reactor safety
transients [@RELAP5VOL1]. It has a great number of two-phase hydraulic
models and heat transfer correlations and has been widely validated with
experiments.

The original version of RELAP does provide a mechanism to handle
calculations coupled to the neutronic code PARCS [@RELAP-PARCS] based in
the PVM libray. However, to be able to perform general coupled
calculations, an extension RELAP5CPL was developed by TECNA S.A. to
exchange information using shared memory objects [@manual-relap-cpl].
Some references that used these extensions include [@relap-fuzzy]
and [@dypra2011]. As RELAP5CPL uses a coupling mechanism compatible with
milonga, it can be used to perform neutronic-thermalhydraulic
calculations.


The problem selected to serve as an example of a coupled calculation is
a simple one-dimensional model both in the thermalhydraulic and in the
neutronic problem. The RELAP input---that is not shown---models a single
pipe subject to a fixed pressure difference. It has some nodes that act
both as a inlet thermalhydraulic zone and as a lower reflector. The
active length consists of twenty volumes that are in contact to a heat
structure that represents the fuel elements. An upper reflector zone is
also provided. As requested in the coupling file---that is shown in the
terminal window---RELAP5CPL writes the twenty values of the fuel
temperature, coolant temperature and coolant density to shared memory
objects. This information is read by milonga and interpolated to
evaluate the cross sections associated to those distributions and to
compute the flux and power distribution, that is written back for
RELAP5CPL to compute again the thermalhydraulics.

The thermalhydraulic model and the neutron cross sections dependance are
made up in order to obtain results that serve to illustrate the effects
of an iterative calculation. The effects of the temperatures and the
xenon on the macroscopic cross sections are exaggerated.

![image](examples//05-coupled/01-relap/relap.pdf)

#### core.mil {#core.mil .unnumbered .unnumbered}

While milonga solves the steady-state neutron diffusion equation, RELAP
is essentially a transient code. In the coupled scheme, RELAP exchange
information once every one hundred transient time steps, that correspond
to a single static step of milonga. The objective of the coupled
calculation is to find a power distribution that generate certain
thermalhydraulic distributions that in turn give rise to the power
distribution found. The details about the thermalhydraulic model are not
important. The coupling file used is shown in the terminal, from which
the general idea of the information exchanged may be grasped. Milonga
reads three vector from shared memory and constructs three continuous
functions containing the coolant temperature, coolant density and fuel
temperature distributions as a function of the axial coordinate $x$.
Conversely, from the continuous `power_density` function it builds a
vector of size twenty that is exported to the shared memory segment by
evaluating the continuous power density at the location of the center of
the RELAP cells measured in milonga's coordinates. Milonga also import a
number of administrative variables written by RELAP, one of which is a
flag to indicate that the calculation has finished. Thus, the number of
static steps is set to a big value to prevent milonga from finishing
before RELAP. The macroscopic nuclear parameters of the reflectors are
constant and uniform, while the cross sections of the active length are
assumed to depend on the fuel burnup distribution, on the fuel
temperature, on the coolant temperature, on the coolant density and on
the xenon distribution in the form

$$\Sigma(x) = \Sigma_0 \big(b(x)\big) + \sum \left. \frac{\partial \Sigma\big(b(x)\big)}{\partial \mathcal{P}} \right|_{\mathcal{P}_0} \cdot \left( \mathcal{P} - \mathcal{P}_0 \right)$$
where $\Sigma_0$ is the cross section evaluated at the nominal
parameters, $b(x)$ is the fuel burnup at position $x$
and $\partial \Sigma/\partial \mathcal{P}$ is the partial derivative of
the cross section with respect to parameter $\mathcal{P}$ evaluated at
the nominal parameter $\mathcal{P}_0$.

The fuel burnup is given by the file `burnup.dat` and the cross sections
and the derivatives are read from `xs.dat`.

The terminal window shows one way of executing this coupled calculation,
namely running RELAP in background by passing an ampersand '&' in the
commandline and then executing milonga afterward. Another way may be to
open two terminals and run each program in one terminal as usual. The
resulting figures show that RELAP start with flat cold temperature
distributions and how they converge to their final values as time passes
by. The exaggerated effects of the thermalhydraulic parameters on the
cross sections, i.e. big derivatives, allow to see how the non-linear
iterative scheme in milonga makes the distributions to oscillate before
converging to the final solution.

\lstset{language=mil,backgroundcolor=\color{mil_fondo},basicstyle=\footnotesize}
    # 1D core coupled with a 1D core model in RELAP
    # a modified version of RELAP called RELAP5CPL is used
    # to couple the codes through shared memory objects and
    # synchronize them using semaphores
    # see references in the documentation for further information

    PROBLEM DIMENSIONS 1 GROUPS 2

    # number of cells in the RELAP nodalization of the core
    relap_number_of_cells = 20

    # this vector contains the coordinates of the centers
    # of the cells that represent the active length in RELAP
    VECTOR relap_cell          20

    # these vectors contain the values that the properties
    # take at the locations given by relap_cell
    # these are to be read from RELAP
    VECTOR coolant_temperature 20
    VECTOR coolant_density     20
    VECTOR fuel_temperature    20

    # construct three continuous functions using the cells
    # positions and the three vectors above that are to be
    # read from RELAP
    FUNCTION Tcool(x) VECTORS relap_cell coolant_temperature INTERPOLATION akima
    FUNCTION Tfuel(x) VECTORS relap_cell fuel_temperature    INTERPOLATION akima
    FUNCTION dcool(x) VECTORS relap_cell coolant_density     INTERPOLATION akima

    # this vector is where the power distribution is to
    # be written by milonga for RELAP to read it
    VECTOR power_distribution  20

    # these variables are used to exchange administrative
    # information with RELAP
    VAR nstsp count iscallrstrec0 problemtype05 problemopt611 restartnum

    # a large number is given here because we expect
    # RELAP5CPL to tell us to stop using variable done
    static_iterations = 1000

    # geometry definition (in cm)
    bot_refl_length = 45
    top_refl_length = 25
    active_length = 530
    x_bare_length = bot_refl_length + active_length + top_refl_length

    # note that the number of cells in the neutronic problem
    # does not need to be the same as in the thermalhydraulic model
    # actually, here they have an offset
    x_cells = 100

    # position of the centers of the active length cells in relap
    # measured in milonga's coordinate system
    relap_cell_1 = 0.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_2 = 1.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_3 = 2.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_4 = 3.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_5 = 4.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_6 = 5.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_7 = 6.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_8 = 7.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_9 = 8.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_10 = 9.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_11 = 10.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_12 = 11.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_13 = 12.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_14 = 13.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_15 = 14.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_16 = 15.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_17 = 16.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_18 = 17.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_19 = 18.5*active_length/relap_number_of_cells + bot_refl_length
    relap_cell_20 = 19.5*active_length/relap_number_of_cells + bot_refl_length

    # the power setpoint for milonga should be in w/cm^2 for
    # for the xenon to be correctly computed
    power = 6500

    # this variable is the actual total power in watts that
    # is needed to convert from power density to the actual
    # power RELAP needs
    total_power = 5e6

    # burnup is given as a function of the axial location measured
    # from the bottom of the active length, to convert it to
    # milonga's coordinate system, the length of the reflector is
    # to be taken into account so b(x) is the burnup to be used
    # care has to be taken only to evaluate b(x) at x located in the
    # active region, otherwise extrapoalated values will be wrong
    FUNCTION b_from_core(x) FILE burnup.dat
    FUNCTION b(x) = b_from_core(x-bot_refl_length)

    # nominal parameters values for the XS tables
    Tfuel0 = 637.0+273.15   # [ K ]
    dcool0 = 801.9          # [ kg / m^3 ]
    Tcool0 = 295.8+273.15   # [ K ]

    # cross sections at central values
    # variable bu means burnup, but it is a dummy variable just
    # to tell milonga these are single-value functions
    # when these funtions are to be evaluated at the MATERIAL
    # keyword, the argument should be b(x)
    FUNCTION D1(bu)    FILE xs.dat COLUMNS 1 2
    FUNCTION D2(bu)    FILE xs.dat COLUMNS 1 3
    FUNCTION abs1(bu)  FILE xs.dat COLUMNS 1 4
    FUNCTION sca12(bu) FILE xs.dat COLUMNS 1 5
    FUNCTION sca21(bu) FILE xs.dat COLUMNS 1 6
    FUNCTION abs2(bu)  FILE xs.dat COLUMNS 1 7
    FUNCTION nuf1(bu)  FILE xs.dat COLUMNS 1 8
    FUNCTION nuf2(bu)  FILE xs.dat COLUMNS 1 9
    FUNCTION Ef1(bu)   FILE xs.dat COLUMNS 1 10
    FUNCTION Ef2(bu)   FILE xs.dat COLUMNS 1 11

    # derivatives with respect to xenon
    FUNCTION dabs1dXe(bu)  FILE xs.dat COLUMNS 1 14
    FUNCTION dsca12dXe(bu) FILE xs.dat COLUMNS 1 15
    FUNCTION dsca21dXe(bu) FILE xs.dat COLUMNS 1 16
    FUNCTION dabs2dXe(bu)  FILE xs.dat COLUMNS 1 17
    FUNCTION dnuf1dXe(bu)  FILE xs.dat COLUMNS 1 18
    FUNCTION dnuf2dXe(bu)  FILE xs.dat COLUMNS 1 19
    FUNCTION dEf1dXe(bu)   FILE xs.dat COLUMNS 1 20
    FUNCTION dEf2dXe(bu)   FILE xs.dat COLUMNS 1 21

    # with respect to fuel temperature
    FUNCTION dabs1dTfuel(bu)  FILE xs.dat COLUMNS 1 24
    FUNCTION dsca12dTfuel(bu) FILE xs.dat COLUMNS 1 25
    FUNCTION dsca21dTfuel(bu) FILE xs.dat COLUMNS 1 26
    FUNCTION dabs2dTfuel(bu)  FILE xs.dat COLUMNS 1 27
    FUNCTION dnuf1dTfuel(bu)  FILE xs.dat COLUMNS 1 28
    FUNCTION dnuf2dTfuel(bu)  FILE xs.dat COLUMNS 1 29
    FUNCTION dEf1dTfuel(bu)   FILE xs.dat COLUMNS 1 30
    FUNCTION dEf2dTfuel(bu)   FILE xs.dat COLUMNS 1 31

    # with respect to coolant temperature
    FUNCTION dabs1dTcool(bu)  FILE xs.dat COLUMNS 1 34
    FUNCTION dsca12dTcool(bu) FILE xs.dat COLUMNS 1 35
    FUNCTION dsca21dTcool(bu) FILE xs.dat COLUMNS 1 36
    FUNCTION dabs2dTcool(bu)  FILE xs.dat COLUMNS 1 37
    FUNCTION dnuf1dTcool(bu)  FILE xs.dat COLUMNS 1 38
    FUNCTION dnuf2dTcool(bu)  FILE xs.dat COLUMNS 1 39
    FUNCTION dEf1dTcool(bu)   FILE xs.dat COLUMNS 1 40
    FUNCTION dEf2dTcool(bu)   FILE xs.dat COLUMNS 1 41

    # with respect to coolant density
    FUNCTION dabs1ddcool(bu)  FILE xs.dat COLUMNS 1 44
    FUNCTION dsca12ddcool(bu) FILE xs.dat COLUMNS 1 45
    FUNCTION dsca21ddcool(bu) FILE xs.dat COLUMNS 1 46
    FUNCTION dabs2ddcool(bu)  FILE xs.dat COLUMNS 1 47
    FUNCTION dnuf1ddcool(bu)  FILE xs.dat COLUMNS 1 48
    FUNCTION dnuf2ddcool(bu)  FILE xs.dat COLUMNS 1 49
    FUNCTION dEf1ddcool(bu)   FILE xs.dat COLUMNS 1 50
    FUNCTION dEf2ddcool(bu)   FILE xs.dat COLUMNS 1 51

    # read RELAP5CPL administrative information including variable
    # done that when set to nonzero makes milonga stop
    PRINT HEADER TEXT "\# waiting for initial semaphore from relap5cpl..." NONEWLINE FLUSH
    IMPORT SHM_OBJECT relap_status {
     SEMAPHORE_WAIT relap_ready
     t dt nstsp count iscallrstrec0 done problemtype05 problemopt611 restartnum
    }
    PRINT HEADER TEXT "got it!"


    # read temperature and density distributions before
    # solving the difussion equation, i.e. before the first
    # ZONE keyword
    IMPORT SHM_OBJECT coolant_temperature coolant_temperature  
    IMPORT SHM_OBJECT coolant_density     coolant_density
    IMPORT SHM_OBJECT fuel_temperature    fuel_temperature

    # give the user some feedback...
    PRINT HEADER TEXT "\# computing initial flux distribution..." NONEWLINE FLUSH

    # core geometry 
    ZONE fuel        MATERIAL fuel
    ZONE refl_bottom MATERIAL refl_bottom X_MAX bot_refl_length
    ZONE refl_top    MATERIAL refl_top    X_MIN x_bare_length-top_refl_length
    PRINT HEADER TEXT "done!"


    # after computing the power distribution (i.e. after the
    # last ZONE) fill in the vector that has to be passed to RELAP
    # relap expects power_distribution to be normalized to one and
    # the total dimensional power in watts as a separate variable
    f = active_length/relap_number_of_cells * 1/power
    power_distribution_1 = power_density(relap_cell_1) * f
    power_distribution_2 = power_density(relap_cell_2) * f
    power_distribution_3 = power_density(relap_cell_3) * f
    power_distribution_4 = power_density(relap_cell_4) * f
    power_distribution_5 = power_density(relap_cell_5) * f
    power_distribution_6 = power_density(relap_cell_6) * f
    power_distribution_7 = power_density(relap_cell_7) * f
    power_distribution_8 = power_density(relap_cell_8) * f
    power_distribution_9 = power_density(relap_cell_9) * f
    power_distribution_10 = power_density(relap_cell_10) * f
    power_distribution_11 = power_density(relap_cell_11) * f
    power_distribution_12 = power_density(relap_cell_12) * f
    power_distribution_13 = power_density(relap_cell_13) * f
    power_distribution_14 = power_density(relap_cell_14) * f
    power_distribution_15 = power_density(relap_cell_15) * f
    power_distribution_16 = power_density(relap_cell_16) * f
    power_distribution_17 = power_density(relap_cell_17) * f
    power_distribution_18 = power_density(relap_cell_18) * f
    power_distribution_19 = power_density(relap_cell_19) * f
    power_distribution_20 = power_density(relap_cell_20) * f

    # write the dimensional power into shared memory
    EXPORT SHM_OBJECT total_power        total_power
    # write the nondimensional power_distribution into shared memory
    # and tell RELAP we are done with our part for this step
    EXPORT SHM_OBJECT power_distribution power_distribution  SEMAPHORE_READY relap_go

    PRINT HEADER TEXT "\# information written and initial semaphore set"

    # print some information to the screen
    PRINT static_step t keff Tcool(relap_cell_20)

    # dump instantaneous distributions to a file
    FILE dist dist.dat STEP
    PRINT_FUNCTION FILE dist power_density xenon Tcool Tfuel dcool SigmaA_1 nuSigmaF_2 flux_1 flux_2

    # be polite, and say goodbye to the user when finished
    PRINT FOOTER TEXT "\# coupled calculation finished, have a nice day!"

    ## material definitions
    # fuel XS are written as a central value plus derivatives times
    # increments with respect to the central values
    MATERIAL fuel {
      D_1         D1(b(x))
      D_2         D2(b(x))
      SigmaA_1    "abs1(b(x)) + dabs1dXe(b(x))*xenon(x) + dabs1dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dabs1dTcool(b(x))*(Tcool(x)-Tcool0) + dabs1ddcool(b(x))*(dcool(x)-dcool0)"
      SigmaA_2    "abs2(b(x)) + dabs2dXe(b(x))*xenon(x) + dabs2dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dabs2dTcool(b(x))*(Tcool(x)-Tcool0) + dabs2ddcool(b(x))*(dcool(x)-dcool0)"

      SigmaS_1.2  "sca12(b(x)) + dsca12dXe(b(x))*xenon(x) + dsca12dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dsca12dTcool(b(x))*(Tcool(x)-Tcool0) + dsca12ddcool(b(x))*(dcool(x)-dcool0)"
      SigmaS_2.1  "sca21(b(x)) + dsca21dXe(b(x))*xenon(x) + dsca21dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dsca21dTcool(b(x))*(Tcool(x)-Tcool0) + dsca21ddcool(b(x))*(dcool(x)-dcool0)"

      nuSigmaF_1  "nuf1(b(x)) + dnuf1dXe(b(x))*xenon(x) + dnuf1dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dnuf1dTcool(b(x))*(Tcool(x)-Tcool0) + dnuf1ddcool(b(x))*(dcool(x)-dcool0)"
      nuSigmaF_2  "nuf2(b(x)) + dnuf2dXe(b(x))*xenon(x) + dnuf2dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dnuf2dTcool(b(x))*(Tcool(x)-Tcool0) + dnuf2ddcool(b(x))*(dcool(x)-dcool0)"

      ESigmaF_1   "1.6e-13*(Ef1(b(x)) + dEf1dXe(b(x))*xenon(x) + dEf1dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dEf1dTcool(b(x))*(Tcool(x)-Tcool0) + dEf1ddcool(b(x))*(dcool(x)-dcool0))"
      ESigmaF_2   "1.6e-13*(Ef2(b(x)) + dEf2dXe(b(x))*xenon(x) + dEf2dTfuel(b(x))*(Tfuel(x)-Tfuel0) + dEf2dTcool(b(x))*(Tcool(x)-Tcool0) + dEf2ddcool(b(x))*(dcool(x)-dcool0))"
    }


    MATERIAL refl_bottom {
      D_1          1.35
      D_2          9.294e-01
      SigmaA_1     7.828e-05
      SigmaA_2     1.463e-02
      SigmaS_1.2   4.368e-05
      SigmaS_2.1   6.507e-04
    }

    MATERIAL refl_top {
      D_1          1.35
      D_2          9.583e-01
      SigmaA_1     1.510e-04
      SigmaA_2     1.234e-02
      SigmaS_1.2   6.891e-05
      SigmaS_2.1   1.059e-03
    }

\lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
    $ cat channel.cpl
    RELAP_TIME_STATUS   relap_status

    RELAP_EXPORT {
     SHARE_NAME coolant_temperature
     SCALAR tempf 87 3:22
     SKIP_STEP 100
    }

    RELAP_EXPORT {
     SHARE_NAME coolant_density
     SCALAR rho    87 3:22
     SKIP_STEP 100
    }

    RELAP_EXPORT {
     SHARE_NAME fuel_temperature
     SCALAR htvatp 1875 1:20
     SEMAPHORE_READY relap_ready
     SKIP_STEP 100
    }

    RELAP_IMPORT {
     SHARE_NAME total_power
     SEMAPHORE_WAIT relap_go
     SCALAR cnvarn 100 0
     SKIP_STEP 100
    }

    RELAP_IMPORT {
     SHARE_NAME power_distribution
     SCALAR cnvsan 201 1
     SCALAR cnvsan 202 1
     SCALAR cnvsan 203 1
     SCALAR cnvsan 204 1
     SCALAR cnvsan 205 1
     SCALAR cnvsan 206 1
     SCALAR cnvsan 207 1
     SCALAR cnvsan 208 1
     SCALAR cnvsan 209 1
     SCALAR cnvsan 210 1
     SCALAR cnvsan 211 1
     SCALAR cnvsan 212 1
     SCALAR cnvsan 213 1
     SCALAR cnvsan 214 1
     SCALAR cnvsan 215 1
     SCALAR cnvsan 216 1
     SCALAR cnvsan 217 1
     SCALAR cnvsan 218 1
     SCALAR cnvsan 219 1
     SCALAR cnvsan 220 1
     SKIP_STEP 100
    }

    RELAP_IMPORT {
     SHARE_NAME relap_done
     SCALAR done 0 0
    }


    $ ./relap5cpl.x -i channel.inp -O channel.out -R channel.rst -a channel.cpl &
    $ milonga core.mil
    ## waiting for initial semaphore from relap5cpl...	relap5cpl: Extension de RELAP5/MOD3.3patch3 para acoplar codigos externos
    relap5cpl: gtheler@tecna.com
    relap5cpl: ultima modificacion  2011-05-02 08:40:14
    relap5cpl: Fecha de compilacion 2011-05-02 08:49:36 linuxgf(jeremy@barnie:x86_64)

    relap5cpl: esperando semaforo "/relap_go"... got it!	
    ## computing initial flux distribution...	ok!
    0======== Execute file name = ./relap5cpl.x                                                                   
                Input file name = channel.inp                                                                             

     Copyright (C) 2001-2006  Information Systems Laboratories, Inc.

     Thermodynamic properties files used by this problem:

     Thermodynamic properties file for d2o      obtained from lfn tpfd2o,
     tpfd2o version 1.0.1, tables of thermodynamic properties of heavy water         
     generated on 09-Jun- 9 at 15:10:20 by stgd2o 1.0 (07/22/91)                     

    0$$$$$$$$ Input processing completed successfully.
      RELAP5/3.3gl          Reactor Loss Of Coolant Analysis Program                                
     Copyright (C) 2001-2006  Information Systems Laboratories, Inc.
    =                                                                                 16-Jul-11    13:13:07      

     cpuT_(s)  probTime  dTime_(s) dTCournt VolCoP  PresCo_MPa VoidCo  emass_kg  VolEr PresEr_MPa  VoidEr    QualaEr   NSteps   Reason
          0.0  0.0000     1.00E-04  0.0     03301    12.14     0.0000   0.00     03301   12.14   0.00    g  0.00    a       0         
    Transient terminated by end of time step cards.
     At time     100.019     seconds; Step      3852
    done!	
    ## information written and initial semaphore set	
    1.000000e+00	0.000000e+00	1.063907e+00	5.510000e+02	
    2.000000e+00	1.397503e+00	8.848904e-01	5.614382e+02	
    3.000000e+00	4.080064e+00	8.939054e-01	5.699261e+02	
    4.000000e+00	6.763396e+00	8.909174e-01	5.748820e+02	
    5.000000e+00	9.425237e+00	8.899390e-01	5.775089e+02	
    6.000000e+00	1.207100e+01	8.891607e-01	5.792419e+02	
    7.000000e+00	1.470793e+01	8.887985e-01	5.801273e+02	
    8.000000e+00	1.733898e+01	8.884757e-01	5.807304e+02	
    9.000000e+00	1.996726e+01	8.883491e-01	5.810124e+02	
    1.000000e+01	2.259325e+01	8.882032e-01	5.812667e+02	
    1.100000e+01	2.521825e+01	8.881619e-01	5.813650e+02	
    1.200000e+01	2.784223e+01	8.880915e-01	5.814805e+02	
    1.300000e+01	3.046586e+01	8.880812e-01	5.815115e+02	
    1.400000e+01	3.308902e+01	8.880452e-01	5.815673e+02	
    1.500000e+01	3.571208e+01	8.880456e-01	5.815744e+02	
    1.600000e+01	3.833490e+01	8.880262e-01	5.816030e+02	
    1.700000e+01	4.095771e+01	8.880296e-01	5.816019e+02	
    1.800000e+01	4.358040e+01	8.880186e-01	5.816174e+02	
    1.900000e+01	4.620310e+01	8.880222e-01	5.816143e+02	
    2.000000e+01	4.882572e+01	8.880157e-01	5.816231e+02	
    2.100000e+01	5.144837e+01	8.880187e-01	5.816200e+02	
    2.200000e+01	5.407098e+01	8.880147e-01	5.816253e+02	
    2.300000e+01	5.669360e+01	8.880170e-01	5.816227e+02	
    2.400000e+01	5.931620e+01	8.880145e-01	5.816260e+02	
    2.500000e+01	6.193881e+01	8.880161e-01	5.816241e+02	
    2.600000e+01	6.456141e+01	8.880145e-01	5.816261e+02	
    2.700000e+01	6.718402e+01	8.880156e-01	5.816248e+02	
    2.800000e+01	6.980661e+01	8.880146e-01	5.816261e+02	
    2.900000e+01	7.242922e+01	8.880154e-01	5.816251e+02	
    3.000000e+01	7.505181e+01	8.880147e-01	5.816260e+02	
    3.100000e+01	7.767441e+01	8.880152e-01	5.816254e+02	
    3.200000e+01	8.029701e+01	8.880147e-01	5.816259e+02	
    3.300000e+01	8.291961e+01	8.880151e-01	5.816255e+02	
    3.400000e+01	8.554221e+01	8.880148e-01	5.816259e+02	
    3.500000e+01	8.816481e+01	8.880150e-01	5.816256e+02	
    3.600000e+01	9.078741e+01	8.880148e-01	5.816258e+02	
    3.700000e+01	9.341001e+01	8.880150e-01	5.816256e+02	
    3.800000e+01	9.603261e+01	8.880149e-01	5.816258e+02	
    3.900000e+01	9.865521e+01	8.880150e-01	5.816256e+02	
    4.000000e+01	1.000190e+02	8.880149e-01	5.816258e+02	
    ## coupled calculation finished, have a nice day!	
    $ gnuplot dists.gnuplot
    $

![image](examples//05-coupled/01-relap/flux.pdf)

![image](examples//05-coupled/01-relap/xenon.pdf)

![image](examples//05-coupled/01-relap/sigmaa1.pdf)

![image](examples//05-coupled/01-relap/nusigmaf2.pdf)

![image](examples//05-coupled/01-relap/power.pdf)

![image](examples//05-coupled/01-relap/coolant.pdf)

![image](examples//05-coupled/01-relap/fuel.pdf)

![image](examples//05-coupled/01-relap/dens.pdf)

\bibliographystyle{unsrt}
Installation and execution
==========================


\hspace{\fill}
\sf 
\footnotesize
Part of the inhumanity of the computer is that, once it is\
competently programmed and working smoothly,\
it is completely honest.


*Isaac Asimov, Change!, 1983*


This chapter gives instructions to perform what should be a one-time
only procedure, namely to compile and install milonga and its required
libraries from scratch. This is the reason why it comes at last: to
avoid appearing in the middle of what should be a periodic-consultation
reference i.e. chapters [2](#cap:equations){reference-type="ref"
reference="cap:equations"}, [3](#cap:input){reference-type="ref"
reference="cap:input"} and [4](#cap:examples){reference-type="ref"
reference="cap:examples"}.

If you are much more impatient than the average user, you can get the
binary versions for your particular architecture
(section [5.1](#sec:veryquick){reference-type="ref"
reference="sec:veryquick"}) instead of trying to compile the source code
(section [5.2](#sec:quick){reference-type="ref" reference="sec:quick"}).
However, it is recommended that you compile the source code to optimize
the execution, especially according to the detailed installation
instructions given in
section [5.3](#sec:detailed-installation){reference-type="ref"
reference="sec:detailed-installation"} to make the most out of milonga.


Milonga is a computer code that instructs a digital computer to perform
certain mathematical operations in order to obtain some results. One
question arises about the type of computers the code should instruct.
Historically, engineering codes were designed to run in mainframes and
supercomputers. As there were a variety of platforms, computer codes
either had to restrict to a single family of processors and operating
systems or to provide means to conditionally include or remove pieces of
source code in order to be compiled with different tools. Nowadays,
desktop---and even laptop---Personal Computers are commonly used as
engineering platforms. Even though there is a wide variety of models,
compatibility and portability between architectures is far more easy
than twenty or thirty years ago, especially if code is developed
according to accepted standards.

Taking into consideration the actual conditions of processors and
operating systems development, the selected platform for developing and
running milonga is a GNU/Linux box running over an Intel-compatible
processor architecture. Besides being free---reason that should be
enough to choose it as a development platform---the GNU operating system
with the Linux kernel[^4] has reached a higher level of maturity and
efficiency running over a wide variety of modern computers than any
other operating system.

Although portability per-se is not into milonga's design basis, it is
expected on the one hand to be able to run the code in a reasonable
spectrum of digital computers and on the other hand to be able to scale
up with both software and hardware developments during a reasonable time
frame, as stated in section [1.2](#sec:basis){reference-type="ref"
reference="sec:basis"}. Thus, even though GNU/Linux is selected as the
development platform, the code is expected to run into other operating
systems and architectures of interest.


Feedback about installation and executing issues---either positive or
negative---is welcome.

Very quick instructions {#sec:veryquick}
-----------------------

1.  Get the binary package for your platform and uncompress it. See

    <http://ib.cnea.gov.ar/~thelerg/wasora/milonga>

    for the list of supported platforms.

    \lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
         $ wget http://ib.cnea.gov.ar/~thelerg/wasora/downloads/milonga-0.1-linux-amd64.tar.gz
         $ tar xvzf milonga-0.1-linux-amd64.tar.gz

2.  Either copy the executable file `milonga` to the directory where
    your input is or vice-versa

    \lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
         $ cd milonga-0.1-linux-amd64
         $ cp examples/test.mil .

3.  Run the program giving the input file as the first argument

    \lstset{language=,backgroundcolor=\color{bash_fondo},basicstyle=\ttfamily\footnotesize}
         $ ./milonga test.mil
         1.092211e+00
         $

A few comments regarding the binaries and its execution:

-   GNU/Linux binaries include PETSc (with all its related libraries it
    needs), SLEPc and GSL as static libraries, but some other common
    libraries such as libglibc, librt and libgfortran are dynamically
    linked. Any modern distribution of GNU/Linux should be able to
    resolve them. If not, please install the associated package.

-   Windows binaries were compiled using Cygwin. They include PETSc,
    SLEPc statically linked in the executable. All the other needed
    libraries are likend dynamically, but the DLLs are distributed
    inside the package under the terms of the GPLv3+.

-   The exact versions of libraries and compilers used to generate each
    executable are stated in each tarball. Again, it is recommended to
    compile the source code using the actual tools installed in the
    target machine.

-   In all the binary versions, PETSc was compiled without MPI support.
    Therefore, the binary may be directly executed by invoking the
    filename from the shell. No MPI wrapper is needed.

-   The execution of milonga in Windows-based architectures is highly
    discouraged! The binary and its related libraries are not designed
    to run natively, and thus their performance is very poor. Microsoft
    Windows is not designed to run engineering codes as milonga. And
    besides, it is not free software. Please try to run milonga in
    GNU/Linux or other unix-based architectures.

Quick instructions {#sec:quick}
------------------

Quick installation instructions for the impatient-but-not-so-much are
given for Debian-based boxes, where milonga was conceived, developed and
fully tested. Nonetheless, the commands given in this section can be
applied to other UNIX variants and even Cygwin distributions where
milonga is known to compile and run, although no comprehensive testing
was made on these platforms.


The list of commands provided in the following sections should do the
trick in a reasonable working GNU/Linux box, i.e. bash, ability to
compile programs and internet access are assumed. If no internet
connection is available (for example because your company keeps blocking
what they think are "security risks" sites or, even worse, do not allow
to connect your GNU/Linux box to the corporate Windows network because
of security issues (sic), of course the downloading stage can be
replaced by a bare copy of the tarballs from a flash drive. In the same
sense, if no Fortran compiler is available, the PETSc compilation
options can be changed to use C-based linear algebra routines. By the
way, it is always handy to have a Fortran compiler ready to go, as you
never know when you may run across a T-Rex---as one that ruled the
continent in page ---with some (of course non-standard) F77 code for you
to run.

Probably some tuning might needed to cope with different software
versions. The script-kiddie approach should be avoided as this
installation procedure may help to get an image of how milonga depends
on these libraries. Commands are given as reference.

### With root access

If you have access to superuser privileges, then some pre-compiled
packages for some of the required libraries can be easily installed
using the distribution's package manager. Example commands are given for
apt-get in Debian-based distributions, but equivalent commands should
also be available in other GNU/Linux distributions.

First make sure your system is able to generate binaries from C sources
(and Fortran for PETSc), has a Python interpreter (also needed to
compile PETSc), manage tarballs and download files using HTTP:

\lstset{language=}
     # apt-get install gcc gfortran make binutils python tar wget 

Then, install the GNU Scientific Libraries (needed by the wasora
framework) and LAPACK (needed by PETSc) development version and its
related dependencies:

     # apt-get install libgsl0-dev liblapack-dev

For some reason, the binary distributions of PETSc and SLEPc provided in
Debian do not work with milonga, so a manual installation is needed. In
any case, compiling both PETSc and SLEPc gives a lot of flexibility in
terms of configuration options that this is a procedure worth learning.
Nevertheless, the following options should be enough:

     $ cd ~
     $ mkdir libs
     $ cd libs
     $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.1-p8.tar.gz
     $ wget http://www.grycap.upv.es/slepc/download/distrib/slepc-3.1-p6.tgz
     $ tar xvzf petsc-lite-3.1-p8.tar.gz
     $ cd petsc-3.1-p8
     $ ./configure --with-shared=1 --with=mpi=0 --with-x=0
     [...]
     $ export PETSC_DIR=$PWD
     $ export PETSC_ARCH=linux-gnu-c-debug
     $ make
     [...]
     $ cd ..
     $ tar xvzf slepc-3.1p-6.tgz
     $ cd slepc-3.1p-6
     $ export SLEPC_DIR=$PWD
     $ ./configure
     [...]
     $ make
     [...]
     $ cd ~

Finally download, compile and install milonga

     $ wget http://ib.cnea.gov.ar/~thelerg/wasora/downloads/milonga-0.1.tar.gz
     $ tar xvzf milonga-0.1.tar.gz
     $ cd milonga-0.1
     $ ./configure
     [...]
     $ make
     [...]
     $ make install
     $ cd ~

These steps should lead to a binary executable of milonga globally
available (in `/usr/local/bin`) for execution:

    $ milonga
    milonga 0.1
    free nuclear reactor core analysis code

    usage: ./milonga input [replacement arguments | PETSc & SLEPc runtime options]
    $ 

### Without root access

If you are planning to run milonga in a host in which you do not have
superuser privileges, then some changes care has to be taken. First,
probably you will not be abe to install GSL libaries as system-wide
accesible, so you will have to compile them in you home directory.
Lucklily, PETSc is able to automatically download and LAPACK withouth
further user interventino. Finally, you will have to tell milonga where
these libraries are located.

     $ cd ~
     $ mkdir libs
     $ cd libs
     $ wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz
     $ wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/fblaslapack-3.1.1.tar.gz
     $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.1-p8.tar.gz
     $ wget http://www.grycap.upv.es/slepc/download/distrib/slepc-3.1-p6.tgz
     $ tar xvzf gsl-1.15.tar.gz
     $ cd gsl-1.15
     $ ./configure
     [...]
     $ make
     [...]
     $ cd ..
     $ tar xvzf petsc-lite-3.1-p8.tar.gz
     $ cd petsc-3.1-p8
     $ ./configure --with-shared=1 --with=mpi=0 --with-x=0 --download-f-blas-lapack
     [...]
     $ export PETSC_DIR=$PWD
     $ export PETSC_ARCH=linux-gnu-c-debug
     $ make
     [...]
     $ cd ..
     $ tar xvzf slepc-3.1p-6.tgz
     $ cd slepc-3.1p-6
     $ export SLEPC_DIR=$PWD
     $ ./configure
     [...]
     $ make
     [...]
     $ cd ~
     $ wget wget http://ib.cnea.gov.ar/~thelerg/wasora/milonga-0.1.tar.gz
     $ tar xvzf milonga-0.1.tar.gz
     $ cd milonga-0.1
     $ ./configure CFLAGS="-I$HOME/libs/gsl-1.15" LDFLAGS="-L$HOME/libs/gsl-1.15/cblas/.libs -L$HOME/libs/gsl-1.15/.libs -Wl,-rpath $HOME/libs/gsl-1.15/cblas/.libs -Wl,-rpath $HOME/libs/gsl-1.15/.libs" 
     [...]
     $ make

You will not be able to install the executable in the system binaries
directory, but you can copy `milonga` to `$HOME/bin` and add this
directoy to the `PATH` environment variable:

     $ mkdir $HOME/bin
     $ cp milonga $HOME/bin
     $ export PATH=$PATH:$HOME/bin
     $ cd
     $ milonga
    milonga 0.1
    free nuclear reactor core analysis code

    usage: ./milonga input [replacement arguments | PETSc & SLEPc runtime options]
     $

Detailed installation instructions {#sec:detailed-installation}
----------------------------------

To make the most out of milonga and to understand how it works, a full
installation of the required libraries is recommended. Users familiar
with GSL, PETSc and SLEPc may want to re-use their installations or to
tune configuration options.

### Obtaining the package

Milonga should be available to download from the author's webpage
located at

<http://ib.cnea.gov.ar/~thelerg/wasora/milonga>

Even though milonga can be redistributed under the terms of the GNU
General Public License[@gpl], the aforementioned site should be the
official source of distribution. Not only should this site contain the
distribution packages milonga version 0.1 but also further news and
related references.

### Required libraries

As discussed in the design basis (), milonga relies on particular
libraries that implement and perform most---if not all---of the
mathematical and numerical operations involved in solving the neutron
diffusion equation. In particular, matrix handling is done by PETSc and
the eigenvalue problem is solved by SLEPc. General mathematical
functions are provided by GSL and CUBATURE is used for multidimensional
adaptive integration.

#### GSL

The GNU Scientific Library is a free libary that provides a great deal
of routines for scientific and engineering computation. At least version
1.14 is needed. It can be obtained from

<http://www.gnu.org/software/gsl/>

To compile milonga, not only the compiled library is needed but also the
headers should be installed where milonga can find them. In most
GNU/Linux distributions, there are packages that provide the needed
functionality. In Debian-based distributions, installing package
`libgsl0-dev` is enough to be able to compile milonga.

     # apt-get install libgsl0-dev

If there is no available package, GSL should be manually installed from
its source code following the standard installation instructions. To
install headers system-wide and make them automatically available to
milonga, you have to be sure to execute `make install` as root after the
compilation:

     $ ./configure
    [...]
     $ make
    [...]
     $ su
    Password: 
     # make install
     # exit
     $

If you do not have root access, either tell GSL's configure script to
install the includes in a directory where you have write permissions and
add it to the `INCLUDE` environment variable or keep the tree where the
GSL source was uncompressed and configure milonga with the appropriate
flags to tell the compiler where the headers are (see
section [5.3.3](#sec:compiling-milonga){reference-type="ref"
reference="sec:compiling-milonga"}). Or, easier, call your system
administrator and kindly ask her to install `libgsl0-dev` from the
distribution repository.


GSL is distributed under the terms of the GNU General Public
License [@gpl].

#### PETSc

PETSc, pronounced PET-see (the S is silent), is a suite of data
structures and routines for the scalable (parallel) solution of
scientific applications modeled by partial differential equations. It
employs the MPI standard for parallelism.


Milonga uses it as an efficient handler of large sparse matrices.
Current version of milonga does not use any of the parallelization
mechanisms provided. At least version 3.1 is needed. Even though most
GNU/Linux distributions do provide a PETSc package, the author was not
able to successfully compile milonga with them. The PETSc source package
distribution has to be downloaded, configured and installed manually in
order to milonga to work.

PETSc's webpage is located at

<http://www.mcs.anl.gov/petsc/petsc-as/>

Detailed installation instructions can be found in the package
documentation [@petsc-user-ref]. A configuration script (written in
Python, so a interpreter has to be installed) detects available
libraries and reads commandline options to generate appropriate
makefiles. Linear algebra libraries and headers LAPACK or BLAS have to
be installed (provided by package `liblapack-dev`, `libblas-dev` or
`libatlas-dev`), although the script can automatically download a free
version of BLAS+LAPACK. It can also detect and use Intel's Math Kernel
Libraries.

It is possible to have a number of PETSc libraries using different
configuration options (using LAPACK or ATLAS, with debug symbols or
optimized, with single or double precision, etc) a target name can be
provided. If no name is provided, default is `linux-gnu-c-debug`.

As milonga 0.1 does not take advantage of PETSc's parallelization
mechanisms, it is best to avoid messing up with MPI, at least for a
first installation.

     $ ./configure --with-mpi=0
    ===============================================================================
                 Configuring PETSc to compile on your system                       
    ===============================================================================
    TESTING: alternateConfigureLibrary from PETSc.packages.petsc4py(config/PETSc/packages/petsc4py.py:70)            Compilers:
      C Compiler:         gcc  -g 
      Fortran Compiler:   gfortran  -g  
    Linkers:
      Static linker:   /usr/bin/ar cr
      Dynamic linker:   /usr/bin/ar  
    X11:
      Includes: 
      Library:  -lX11
    BLAS/LAPACK: -llapack -lblas
    PETSc:
      PETSC_ARCH: linux-gnu-c-debug
      PETSC_DIR: /home/jeremy/libs/petsc-3.1-p8
      Clanguage: C
      Scalar type: real
      Precision: double
      Memory alignment: 16
      shared libraries: disabled
      dynamic libraries: disabled
    xxx=========================================================================xxx
       Configure stage complete. Now build PETSc libraries with:
       make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-c-debug all
    xxx=========================================================================xxx
     $ make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-c-debug all
    [...]
    libfast in: /home/jeremy/libs/petsc-3.1-p8/tutorials/multiphysics
    Completed building libraries
    =========================================
    Now to check if the libraries are working do:
    make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-c-debug test
    =========================================
     $ make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-c-debug test
    Running test examples to verify correct installation
    C/C++ example src/snes/examples/tutorials/ex19 run successfully with 1 MPI process
    Fortran example src/snes/examples/tutorials/ex5f run successfully with 1 MPI process
    Completed test examples
     $

If `configure` complains about missing libraries, further reductions are
possible. Usually, graphical X routines are not needed, so `configure`
can be told to automatically download algebra libraries. Also, if no
debugging is needed, an optimized library can be compiled.

     $ ./configure --with-mpi=0 --with-x=0 --download-f-blas-lapack=1 --with-debugging=0
    [...]
     $

Additionally, shared and/or dynamic versions can be generated by using
`–with-shared` and `–with-dynamic`. This way, the milonga executable
will be very much smaller than with the default static PETSc
configuration. If no Fortran compiler is available,
`–download-c-blas-lapack` should be given to instruct the script to
download, compile and use a BLAS+LAPACK version converted from Fortran
to C using `f2c`.


To compile SLEPc and milonga, two environment variables have to be set.
The first one is `PETSC_DIR` and should contain the directory where
PETSc was uncompressed. It can be easily set by executing

     $ export PETSC_DIR=$PWD
     $

from the PETSc directory. The other one is called `PETSC_ARCH` and is
the selected target. In GNU/Linux, default is `linux-gnu-c-debug` as
reported by the configure output.

     $ export PETSC_ARCH=linux-gnu-c-debug
     $

Other targets can be generated by passing the `–with-petsc-arch`
argument to `configure`. For example, to compile an optimized shared
version of PETSc using a C-based LAPACK library, configure with

     $ ./configure --with-petsc-arch=linux-gnu-opt --with-debugging=0 --with-mpi=0 --with-x=0 --with-shared=1 --with-fc=0 --download-c-blas-lapack=1
    ===============================================================================
                 Configuring PETSc to compile on your system                       
    ===============================================================================
    [...]
    xxx=========================================================================xxx
       Configure stage complete. Now build PETSc libraries with:
       make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-opt all
    xxx=========================================================================xxx
     $ make PETSC_DIR=/home/jeremy/libs/petsc-3.1-p8 PETSC_ARCH=linux-gnu-opt all
    [...]
     $

and then set `PETSC_ARCH` to the appropriate value

     $ export PETSC_ARCH=linux-gnu-opt
     $


PETSs is released under a free license developed by the University of
Chicago [@petsc-license].

#### SLEPc

SLEPc is a software library for the solution of large scale sparse
eigenvalue problems on parallel computers. It is an extension of PETSc
and can be used for either standard or generalized eigenproblems, with
real or complex arithmetic. it can also be used for computing a partial
SVD of a large, sparse, rectangular matrix, and to solve quadratic
eigenvalue problems.

It is used by milonga to solve the generalized eigenvalue problem
associated to the multigroup steady-state neutron diffusion problem in a
discretized spatial domain. At least version 3.1 is needed. As with
PETSc, packages provided by GNU/Linux distributions do not work with
milonga. SLEPc's webpage is

<http://http://www.grycap.upv.es/slepc/>

For complete instructions on the installation of SLEPc refer to its
documentation[@slepc-users-manual]. It also provides a configuration
script written in Python. Luckily, it reads the configuration options
passed to PETSc so normally no further arguments are needed. However, in
order to be able to read PETSc's configuration, both
variables `PETSC_DIR` and `PETSC_ARCH` have to be defined to reflect the
installation directory and the target architecture. For example, if
PETSc was configured with the default architecture and the steps in the
previous section were followed, the execution of SLEPC's `configure`
should be simple

     $ ./configure
    Checking environment...
    Checking PETSc installation...
    Checking LAPACK library...

    ================================================================================
    SLEPc Configuration
    ================================================================================

    SLEPc source directory:
     /home/jeremy/libs/slepc-3.1-p6
    SLEPc install directory:
     /home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug
    PETSc directory:
     /home/jeremy/libs/petsc-3.1-p8
    Architecture "linux-gnu-c-debug" with double precision real numbers
     $

Before compiling SLEPc, the `SLEPC_DIR` variable has to be set to
reflect SLEPc's directory. It can be easily filled with the current
directory by using the `PWD` variable. A `make` command will compile the
library.

     $ export SLEPC_DIR=$PWD
     $ make
    ==========================================
    On Fri Jul 15 19:48:34 ART 2011 on tom
    Machine characteristics: Linux tom 2.6.32-5-amd64 #1 SMP Mon Mar 7 21:35:22 UTC 2011 x86_64 GNU/Linux
    -----------------------------------------
    Using SLEPc directory: /home/jeremy/libs/slepc-3.1-p6
    Using PETSc directory: /home/jeremy/libs/petsc-3.1-p8
    Using PETSc arch: linux-gnu-c-debug
    -----------------------------------------
    SLEPC_VERSION_RELEASE    1
    [...]
    libfast in: /home/jeremy/libs/slepc-3.1-p6/include/private
    libfast in: /home/jeremy/libs/slepc-3.1-p6/docs
    /usr/bin/ranlib /home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/lib/*.a
    Completed building SLEPc libraries
    =========================================
    making shared libraries in /home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/lib
    building libslepc.so
    =========================================
    Now to check if the libraries are working do: make test
    =========================================
     $ make test
    Running test examples to verify correct installation
    C/C++ example src/examples/ex1 run successfully with 1 MPI process
    Fortran example src/examples/ex1f run successfully with 1 MPI process
    Completed test examples
     $ 

SLEPc is released under the terms of the Lesser General Public
License [@lgpl].

#### cubature

Cubature is a simple C subroutine for adaptive multidimensional
integration of vector-valued integrands over hypercubes. Of course it
can handle scalar integrands as the special cases of one-dimensional
vectors. It can be obtained from

<http://ab-initio.mit.edu/wiki/index.php/Cubature/>

It is used by milonga to compute cell average cross sections and to
compute the surface integrals of the diffusion coefficients. As cubature
is a small routine, it is incorporated into milonga as another source
file so actually no installation is needed.


Cubature is released under the terms of the GNU General Public License
v2 or later [@gpl].

### Compiling milonga {#sec:compiling-milonga}

Lastly, it is time to compile milonga. If GSL is installed in default
the location and the PETSc and SLEPc environment variables are set to
their proper values, then the standard `./configure` and `make` steps
should do the trick:

     $ ./configure
    checking for a BSD-compatible install... /usr/bin/install -c
    checking whether build environment is sane... yes
    checking for a thread-safe mkdir -p... /bin/mkdir -p
    checking for gawk... gawk
    checking whether make sets $(MAKE)... yes
    checking for PETSC_DIR environment variable... /home/jeremy/libs/petsc-3.1-p8
    checking for PETSC_ARCH environment variable... linux-gnu-c-debug
    checking for SLEPC_DIR environment variable... /home/jeremy/libs/slepc-3.1-p6
    checking for gcc... gcc
    [...]
    checking for strdup... yes
    checking for strstr... yes
    configure: creating ./config.status
    config.status: creating Makefile
    config.status: creating src/config.h
    config.status: src/config.h is unchanged
    config.status: executing depfiles commands
     $ make
    gcc -DHAVE_CONFIG_H -I. -I./src     -g -O2 -I/home/jeremy/libs/petsc-3.1-p8/linux-gnu-c-debug/include -I/home/jeremy/libs/petsc-3.1-p8/include -I/home/jeremy/libs/petsc-3.1-p8/include/mpiuni   -I/home/jeremy/libs/slepc-3.1-p6 -I/home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/include -I/home/jeremy/libs/slepc-3.1-p6/include -MT assignment.o -MD -MP -MF .deps/assignment.Tpo -c -o assignment.o `test -f 'src/assignment.c' || echo './'`src/assignment.c
    mv -f .deps/assignment.Tpo .deps/assignment.Po
    [...]
    mv -f .deps/hello.Tpo .deps/hello.Po
    gcc  -g -O2 -I/home/jeremy/libs/petsc-3.1-p8/linux-gnu-c-debug/include -I/home/jeremy/libs/petsc-3.1-p8/include -I/home/jeremy/libs/petsc-3.1-p8/include/mpiuni   -I/home/jeremy/libs/slepc-3.1-p6 -I/home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/include -I/home/jeremy/libs/slepc-3.1-p6/include   -o milonga assignment.o builtinfunctionals.o builtinfunctions.o call.o cleanup.o commoninit.o commonparser.o error.o function.o getptr.o handler.o history.o instruction.o io.o line.o parametric.o parseaux.o print.o realtime.o shmem.o uservars.o allocate.o boundary.o cubature.o debug.o eigen.o geometry.o init.o matrices.o milonga.o outputs.o parser.o petschandler.o power.o printresult.o result.o steady.o step.o version.o xenon.o xs.o hello.o  -lgsl -lgslcblas -lm -lpthread -lrt  -Wl,-rpath,/home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/lib -L/home/jeremy/libs/slepc-3.1-p6/linux-gnu-c-debug/lib -lslepc      -Wl,-rpath,/home/jeremy/libs/petsc-3.1-p8/linux-gnu-c-debug/lib -L/home/jeremy/libs/petsc-3.1-p8/linux-gnu-c-debug/lib -lpetsc   -llapack -lblas -lm -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.5 -Wl,-rpath-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -lgcc_s -lgfortran -lm -lm -ldl -lgcc_s -ldl 
     $

If variables `PETSC_DIR`, `PETSC_ARCH` or `SLEPC_DIR` are not set,
configuration script will complain. However, if they are set but their
content is not correct, `configure` will run but `make` will fail.

If GSL was installed either from the distribution's packages or compiled
from the source and then installed as root, the header should be
globally available and milonga's `configure` should be able to find
them. If, however, they are installed somewhere else---probably because
no root access is available---then you have to tell milonga where they
are. For example, if you uncompressed the source distribution into
directory ` /libs/gsl-1.15` and performed the `./configure` and `make`
commmands there---without the `make install` step---then you have to
pass to milonga's configuration script the following arguments:

     $ ./configure CFLAGS="-I$HOME/libs/gsl-1.15" LDFLAGS="-L$HOME/libs/gsl-1.15/cblas/.libs -L$HOME/libs/gsl-1.15/.libs -Wl,-rpath $HOME/libs/gsl-1.15/cblas/.libs -Wl,-rpath $HOME/libs/gsl-1.15/.libs" 
    [...]
     $ make
    [...]
     $

If you want to have milonga as a globally-available command---and you
have superuser access of course---you can have it installed in
`/usr/local/bin` by issuing

     # make install
     #

Milonga can now be called from any directory now.

Execution
---------

After a successful compilation---or by directly downloading a binary
distribution---you will end up with an executable called milonga.
Scientific codes linked against PETSc usually make use its
parallelization capabilities and, as such, they have to be run using an
MPI wrapper application. However, current version of milonga does not
support parallelization, and thus it can be executed as any other
regular program.

Milonga needs as its first argument the path of the input file
containing the problem it has to solve. If milonga is system-wide
available, you can go to the directory where the input is and just call
milonga followed by the input name

     $ cd examples
     $ milonga test.mil
    1.092211e+00
     $

If milonga is not installed system-wide, either you will have to copy
the executable where the input is, provide the full path to the
executable or provide the full path to the input. Keep in mind that
using the last choice that files that are included using relative paths
will not be found.

For example, the examples of chapter four can be run by using the `runx`
script provided in the `examples` subdirectory of milonga's
distribution:

     $ cd examples
     $ cd 01-analytical
     $ ./runx
    01-bare_slab
    milonga comparisson.mil
    analytical keff =       1.09220381
    numerical  keff =       1.09221091
         difference =       7.101524e-06
    milonga flux.mil
    5.000000e-01    2.466006e-02    2.467300e-02
    1.500000e+00    7.398017e-02    7.399464e-02
    2.500000e+00    1.232273e-01    1.232433e-01
    [...]
    06-two-zone_slab
    milonga twozone.mil
    ## copy and paste the following lines into gnuplot to obtain
    ## the continuous flux distribution as a function of x
    numerical_keff  =       1.14799085
    analytical_keff =       1.14779361
    phi1(x) =  1.316009e-02  * sinh( 1.121395e-01  * x ) 
    phi2(x) =  1.316009e-02  * sinh( 5.606974e+00  ) / sin( 2.261471e+00 ) * sin( 4.522943e-02 *( 1.000000e+02 -x ) ) 
    phi(x) = ( x < 5.000000e+01 ) ? phi1(x) : phi2(x) 


    good! no errors found!
     $

Syntax highlight
----------------

Input files are used to instruct a digital computer to perform a certain
task. As discussed in the design basis in
chapter [1](#cap:introduction){reference-type="ref"
reference="cap:introduction"}, they should in some sense resemble the
approach taken by high-level source code compilers. This way, milonga's
input files should be friendly human-readable chunks of text instead of
a bunch of numbers representing perforated cards.

As with source code, to render input files even more human-friendly, it
is very handy to have a system of syntax highlight when working with
milonga's input files. The distribution includes a syntax file
definition for the text editor kate.

![[\[fig:kate\]]{#fig:kate label="fig:kate"} Syntax highlight in editor
kate for milonga input files.](kate)


To have kate understand milonga's input file, copy `milonga.xml`
(located in subdirectory `doc`) to katepart's syntax directory. You may
have to create the directory if it does not exist:

     $ mkdir -p $HOME/.kde/share/apps/katepart/syntax
     $ cp doc/milonga.xml $HOME/.kde/share/apps/katepart/syntax
     $

This way, any file with extension `mil` opened in kate or in any other
editor using kate's framework (such as kwrite) will be highlighted as
shown in figure [\[fig:kate\]](#fig:kate){reference-type="ref"
reference="fig:kate"}.

\bibliographystyle{unsrt}
\appendix
\backmatter
\printindex
\newpage
\thispagestyle{empty}
 


\
\

\
\

\vspace{1cm}
<http://ib.cnea.gov.ar/~thelerg/wasora/milonga>

\vspace{0.5cm}
![image](sombrero){width="2cm"}

[^1]: As of 2011.

[^2]: For the record, this sentence was first written in 2011.

[^3]: Although the current version does not utilize all the potential
    this library provides, especially the parallelization capabiltities.
    Future versions will be fully optimized with scalability in mind.

[^4]: Strictly speaking, the official Linux kernel release contains some
    portions of code that are not free in the GNU sense. However, it is
    possible to obtain 100% free versions of the kernel to build free
    GNU/Linux distributions such as Debian.
