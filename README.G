/*
                    THE *.G PROCEDURES FROM S. GREENLAND
                             15 Oct. 2007

      This is a set of  GAUSS procedures for statisical
    analysis and stat programming. The procedures range from one-line
    language extensions (such as MAX.G, MIN.G, & SUMR.G) that simplify
    coding, to full programs for nonlinear regression models, including
    a few (such as stereotype logistic regression, LRORDST.G; logistic
    regression for two-stage and nested studies, LR2STAGE.G and LRPSEUDO.G;
    and penalized-likelihood regression, EXPREGP.G, EXPRISKP.G, LOGREGP.G,
    LRPM.G, LRPOLYP.G, and LRPMP.G) that are currently not available in
    most standard packages.
  
      EACH PROCEDURE FILE BEGINS WITH A PROCEDURE DESCRIPTION. READ IT
    TO FIND OUT HOW TO USE THE PROCEDURE. Read the procedure code to see
    how it works.

       This procedure set is experimental and is under continuous
    development. Copies or updates may be obtained by sending an e-mail 
   with your title & address, to Sander Greenland:  lesdomes@ucla.edu
    Incomplete requests will be ignored. If you use these procedures
    for more than a year, an update is strongly recommended.

       If you use or examine them enough, you will undoubtedly discover
    bugs, stat errors, and programming inefficiencies in the procedures.
    I would greatly appreciate being notified of anything you find, and
    would welcome suggestions for improving the code or set. In coding
    these programs, I have tried to strike a balance between transparency
    and efficiency, but in some cases may be far from both. Because the
    regression procedures are partly for instructional purposes, however,
    I have used only basic nonlinear algorithms.

    REQUIREMENTS:
       You must have a licensed copy of GAUSS 3.2 or higher and run the
    procedures from a GAUSS program or the GAUSS interpreter. This means that
    you must know how to read and manipulate your data and call procedures in
    GAUSS.

    INSTALLATION:
       Currently, you have at least two options for mounting the set:
    Option 1: Create a subdirectory of your GAUSS directory called
      (for example) PRG. Copy the entire set into that subdirectory. Then
      use a text editor to modify the GAUSS configuration program
      GAUSSI.CFG in the GAUSS directory as follows:
      Change the line
      src_path = $(GAUSSDIR)\src
      to
      src_path = $(GAUSSDIR)\src;$(GAUSSDIR)\prg
      YOU MUST DO THIS EVERY TIME YOU INSTALL GAUSS OR A GAUSS UPDATE.
    Option 2: Load the set into the SRC subdirectory.

    BUG IN GAUSS COMPILER BEFORE VERSION 4:
      There is a bug in earlier Gauss compilers which causes compilation to
    terminate with a "Wrong number of returns" error when calling certain
    *.G procedures under certain conditions, especially DIM.G and INVNORM.G.
    If this should happen, simply append the cited procedure to the end of
    the calling program. 

    WARNINGS:
       Although the procedures are simple ASCII text files, practice
    safe file transfer protocol: virus scan the disk before installing.

    DISCLAIMER:
       This procedure set is supplied with absolutely no warranties
    expressed or implied. Sander Greenland accepts no responsibility
    for any damages arising from its use.
*/