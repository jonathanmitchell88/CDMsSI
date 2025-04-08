####################
# Article: Convergence-divergence models: Generalizations of phylogenetic trees modeling gene flow over time
# Journal: For peer review
# Authors: Jonathan D. Mitchell (1,2*) and Barbara R. Holland (1,2)
# Affiliations: 1) School of Natural Sciences (Mathematics), University of Tasmania, Hobart, TAS, Australia, 2) ARC Centre of Excellence for Plant Success in Nature and Agriculture, University of Tasmania, Hobart, TAS, Australia.
# Corresponding author e-mail: jonathanmitchell88@gmail.com
####################

-- CDM 5
-- We want to find the set of equations that can be used to solve for the parameters of CDM 5.

-- r_1=r0011
-- r_2=r0101
-- r_3=r0110
-- r_4=r0111
-- r_5=r1001
-- r_6=r1010
-- r_7=r1011
-- r_8=r1100
-- r_9=r1101
-- r_10=r1110
-- r_11=r1111
-- r_12=delta

-- To solve for the x parameters we want a set of rs with as large a cardinality as possible such that no rs are functions of other rs in the set. q1111 depends on both r_11 and r_12, with no other qs depending on r_11 or r_12. Thus, at most only one of r_11 and r_12 can be included in the set. We find the largest cardinality set that includes at most one of r_11 and r_12.

S=QQ
R=S[x_1..x_11,r_1..r_12,MonomialOrder=>{GRevLex=>11,GRevLex=>12}]
I1=ideal(r_1-x_4*x_5*x_6*x_7*x_9*x_10,r_2-x_10*x_11*(1-x_9*(1-x_2*x_3*x_4*x_6*x_8)),r_3-x_7*x_8*x_9*x_11*(1-x_6*(1-x_2*x_3*x_5)),r_4-x_2*x_3*x_4*x_5*x_6*x_7*x_8*x_9*x_10*x_11,r_5-x_1*x_2*x_4*x_9*x_10,r_6-x_1*x_2*x_5*x_6*x_7,r_7-x_1*x_2*x_4*x_5*x_6*x_7*x_9*x_10,r_8-x_1*x_3*x_6*x_8*x_9*x_11,r_9-x_1*x_2*x_3*x_4*x_6*x_8*x_9*x_10*x_11,r_10-x_1*x_2*x_3*x_5*x_6*x_7*x_8*x_9*x_11,r_11-x_1*x_7*x_10*x_11*(x_4*x_8*x_9*(x_2*(1-x_6)+x_3*x_5*x_6)+x_2*x_5*x_6*(1-x_9)),r_12-x_1*x_2*x_3*x_4*x_5*x_6*x_7*x_8*x_9*x_10*x_11)
I2=ideal(gens gb I1)

-- Eliminate x_1..x_11. Since we can include at most one of r_11 and r_12, first try eliminating r_11.

I3=eliminate(I2,{x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10,x_11,r_11})

-- I3 has no generators that are functions of r_12. Thus, r_12 is in the set. (We have not checked whether r_11 could be included in the set without r_12, but this is not necessary.)

-- We have eliminated 12 variables. Thus, the cardinality of the set is:

dim I3-12

-- gamma is another variable to include in the set. Thus, the cardinality of the set including gamma is:

dim I3-12+1

-- This is the set {r_1,r_2,r_3,r_4,r_5,r_6,r_8,r_9,r_12,gamma}.
