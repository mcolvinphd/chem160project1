from random import choice
npart=300
side=35  #Should be an odd number
time=0
steps = [(1,0),(-1,0),(0,1),(0,-1)]
grid=[[0 for x in range(side)] for y in range(side)]
for ipart in range(npart):
    # Start particle at center
    x,y = side//2,side//2
    counter=0
    # perform the random walk until particle departs
    while 1:
        counter+=1
        grid[x][y]=0   #Remove particle from current spot
        # Randomly move particle
        sx,sy = choice(steps)
        x += sx
        y += sy
        if x<0 or y<0 or x==side or y==side:
            time+=counter
            break
        grid[x][y]=1   #Put particle in new location
avetime=time/npart
print("<t>=%5.2f  <t>/r2=%5.2f"%(avetime, avetime/(side**2)))

