#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#%%
# An "interface" to matplotlib.axes.Axes.hist() method
d = np.random.laplace(loc=15, scale=3, size=500)
n, bins, patches = plt.hist(x=d, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
x = np.random.random_integers(1, 100, 5)
plt.hist(x, bins=20)
plt.ylabel('No of times')
plt.show()

# %%
import matplotlib.pyplot as plt
 
x = [1,1,2,3,3,5,7,8,9,10,
     10,11,11,13,13,15,16,17,18,18,
     18,19,20,21,21,23,24,24,25,25,
     25,25,26,26,26,27,27,27,27,27,
     29,30,30,31,33,34,34,34,35,36,
     36,37,37,38,38,39,40,41,41,42,
     43,44,45,45,46,47,48,48,49,50,
     51,52,53,54,55,55,56,57,58,60,
     61,63,64,65,66,68,70,71,72,74,
     75,77,81,83,84,87,89,90,90,91
     ]
y=[1,3,5,15,16,22,22,22,22,22,33]
plt.style.use('ggplot')

#plt.hist(y, bins=[0,10,20,30,40,50,60,70,80,90,99])

n, bins, patches = plt.hist(x=y, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
#plt.text(10, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


# %%

def plotHistogramOfTimeDifferences(filepath, timestampPosition, maxXLabel, binSize):
    timeDifferences=[]
    previousTimestamp=0

    with open(filepath, "r") as f:
        mylist = f.read().splitlines()
        #print('mylist:::')
        #print(mylist)
        for lines in mylist:
            values = lines.split(" ")
            timestamp=int(values[timestampPosition])
            if(previousTimestamp!=0):
                timeDifferences.append(previousTimestamp-timestamp)
            
            previousTimestamp=timestamp


            #print('v:')
            #print(values[timestampPosition])
            #print('')

    #print('')
    print('timeDifferences:')
    print(timeDifferences)

    maxTimeDifferences=max(timeDifferences)
    print(maxTimeDifferences)
    print(type(3))
    print(type(maxTimeDifferences))

    myBins=range(0,min(maxXLabel,maxTimeDifferences)+2,binSize)
    
    #n, bins, patches = plt.hist(x=np.clip(timeDifferences,0,myBins[-1]), bins=myBins, color='#0504aa', alpha=0.7, rwidth=0.85)
    n, bins, patches = plt.hist(x=timeDifferences, bins=myBins, color='#0504aa', alpha=0.7, rwidth=0.85)
    print('n')
    print(n)
    #print(bins)
    #print(patches)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('time difference in seconds')
    plt.ylabel('Frequency')
    plt.title('Histogram of time differences between gps points')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


    
    xlabels = bins[0:].astype(str)
    print('test:')
    print(xlabels)
    xlabels[-1] += '+'
    print('test')
    print(type(xlabels))

    myXTicks=np.array(myBins)+(binSize/2)
    #myXTicks=np.array(myBins)
    print('myXTicks')
    print(myXTicks)
    print('xlabel')
    print(xlabels)
    print('myBins')
    print(myBins)
    plt.xticks(myXTicks,xlabels)
    

    #N_labels = len(xlabels)
    #plt.xlim([0, 325])
    #plt.xticks(25 * np.arange(N_labels) + 12.5)
    #ax.set_xticklabels(xlabels)
    
    



path='/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy_verysmall_closestIsBest_1stSolution.txt'
path2='/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy_small.txt'
path3='/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy_small_verysmall_3.txt'


plotHistogramOfTimeDifferences(path3, 3,15,1)



# %%
x=range(1,10,2)
for x in x:
    print(x)


# %%
