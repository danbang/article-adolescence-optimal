# article-adolescence-optimal

This repository contains anonymised behavioural data and code supporting the following paper:

Haller*, Bang*, Bahrami &amp; Lau, Group-decision making is optimal in adolescence. Under review.

FigureX.m files will generate the specified plots and analyses from the paper.

The "Helpers" folder should be downloaded; it contains a custom script for computing sensitivity.

**The data file contains the following information**

*groupID*: group identifier

*groupAGE*: group mean age

*condition*: age group (1: younger adolescents; 2: older adolescents; 3: adults)

*sbjID*: subject identifier

*sbjNUM*: subject within dyad identifier (1: keyboard in session 1 and mouse in session 2; 2: mouse in session 1 and keyboard in session 2)

*sbjYEARS*: subject age in years

*sbjMONTHS*: subject age in months

*trial*: trial identifier

*session*: session identifier (session 1: trials 1-128; session 2: trials 129-256)

*stimInterval*: target interval

*stimContrast*: target contrast

*stimDelta*: target contrast * interval (negative: interval 1; positive: interval 2)

*sbjChoice*: subject choice 

*sbjAcc*: subject accuracy

*sbjRT*: subject reaction time (only valid when subject is using keyboard - see *sbjNUM*)

*sbjArbi*: subject indicating joint decision if disagreement

*othChoice*: partner choice

*othAcc*: partner accuracy

*othRT*: partner reaction time (only valid when subject is using keyboard - see *sbjNUM*)

*othArbi*: partner indicating joint decision if disagreement

*dyaChoice*: dyad choice

*dyaAcc*: dyad accuracy

*dyaRT*: dyad deliberation time if applicable

*disagree*: disagreement (0: no; 1: yes)

This material is being released with a permissive open-source license. If you make use of the data or ode, we would appreciate that you cite the paper.
