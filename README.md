# Fungal_Spillover_Code
Fungal spillover influences plant competitive ability and affects species’ coexistence

Summary: Tropical forests are the most biodiverse ecosystems in the world and are critical to global ecosystem functions. Anthropogenic activities have increased forest fragmentation, which could introduce new biotic interactions into tropical forests via spillover, yet little is known about how spillover influences tropical tree diversity. Existing spillover theory has mainly focused on transmission dynamics rather than plant species coexistence. We developed a theoretical spatial model consisting of one predator (fungi) interacting with two prey (plant) populations and introduced a new fungal pathogen via diffusion to mimic spillover. We aimed to understand how spillover pathogen and resident pathogen host-specificity, spillover intensity, and plant competitive abilities alter plant coexistence outcomes. Results show that spillover, contingent on plant species’ competitive ability, can either promote or disrupt coexistence, regardless of the host specificity of the resident pathogen. While generalist spillover disrupts plant coexistence at the forest edge, specialized spillover acts as a stabilizing mechanism when targeting the strong plant competitor, eventually reversing plant apparent competitive abilities. These effects are heightened with wind. Our model provides insight on the influences of fungal pathogen spillover on plant community composition and diversity when the effects of fungal spillover in field are unknown.

Data
turb.csv a vector of random turbulence values used to simulate turbulent winds.

Code
Code_final.R is the script associated with the manuscript which uses the software R

Packages needed deSolve, ReacTran, abc, abcrf, BioManager, ggplot2, dplyer, lime, keras, tensorflow, EBImage, tfdatasets, RColorBrewer, lattice, vegan

Brief description of the script
To create a randomized vector to simulate turbulence that will be used in a later model. #lines 18-25

Model0: baseline model with two plant species interacting with a resident pathogen in the absence of spillover. #lines 31-38

Model1: Testing the effects of introduced pathogen. Two plant species interacting with both the resident and introduced pathogen #lines 40-56

Model2: Testing the effects of introduced pathogen with the addition of wind. Two plant species interacting with both the resident and introduced pathogen #lines 58-&4

Model 3: Testing the effects of introduced pathogen with simulated turbulence, using the vector created in lines 18-25. Two plant species interacting with both the resident and introduced pathogen #lines 76-92

##To change wind speed (velocity) 
Wind & no velocity: v=0 #line 52
Wind & velocity: v=0.25 #line 70
Wind & turbulence: v=turb #line 88; (turb= turbulence vector #line 85)

Investigation of the effects of intraspecific competition on the coexistence outcome of model0. #lines 94-135

Extracting the specific population dynamics #lines 137-144

Summarizing and plotting the outcome of intraspecific competition on long-term coexistence: #lines 146-164 P1<=0.00001 P2>0.00001 Y1>0, P2 outcompetes P1. P1<=0.00001 P2>0.00001 Y1>0, P1 outcompetes P2. P1>0.00001 P2>0.00001 Y1>0, Coexistence. P1<=0.00001 P2<=0.00001 Y1<0, System crash.

Repeating for diffusion model 1 #lines 191-278

Setting up the grid space #lines 194-195

Summarizing and plotting the outcome of intraspecific competition on long-term coexistence: #lines 253-277
P1<=0.00001 P2>0.00001 Y1>0, P2 outcompetes P1. P1<=0.00001 P2>0.00001 Y1>0, P1 outcompetes P2. P1>0.00001 P2>0.00001 Y1>0, Coexistence. P1<=0.00001 P2<=0.00001 Y1<0, System crash.

Repeat for model2  #lines 279-355

Finally repeat for model3 #lines 356-432

Finally comment in line 433 explain how everything was repeated to vary between specialist or generalist spillover and resident pathogens.

Performing 10,000 simulations by randomly drawing each parameters. The simulations are then analyzed using redundancy analyses and abc random forests. # lines 436-881
