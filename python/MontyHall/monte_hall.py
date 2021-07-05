import numpy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

import sys
if len(sys.argv)<2: sys.exit("run with python monte_hall.py swap/stick [Ngames]")
prize_boxes = ["A","B","C"]

strategy = sys.argv[1]
if strategy not in ["swap","stick"]: sys.exit("choose one of 'swap', 'stick'")
if len(sys.argv) > 2 : ngames = int(sys.argv[2])
else: ngames = 1

games = []
for n in range(ngames) :

   prize = numpy.random.choice(prize_boxes)

   # player chooses a box
   player_choice = numpy.random.choice(prize_boxes)

   # host opens a box from remaning boxes.
   remaining_choices = [c for c in prize_boxes if c not in [prize,player_choice]]
   host_choice = numpy.random.choice(remaining_choices)
   remaining_box = [c for c in prize_boxes if c not in [host_choice,player_choice]]

   # now player can swap or stick
   if strategy == "swap":
    original_choice = player_choice
    player_choice = remaining_box[0]
    if ngames <=1 : print("Player chose box %s originally, host opened box %s, player swapped to box %s"%(original_choice,host_choice,player_choice))
   else:
    if ngames <=1 : print("Player chose box %s, host opened box %s"%(player_choice,host_choice))
   if ngames <=1 :
     print(" .... prize was in box %s"%(prize))
     if player_choice==prize: print(" -----> WIN!")
     else: print(" -----> LOSE :(")

   if player_choice==prize: games.append("win")
   else : games.append("lose")

if ngames>1:
  plt.hist(games)
  plt.show()
