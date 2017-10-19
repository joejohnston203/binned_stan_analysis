def time_profile_365_bins(x):
    if(0.0<=x<10.0/365.0):
        return 1
    elif(x<20.0/365.0):
        return 2
    elif(x<90.0/365.0):
        return 1
    elif(x<100.0/365.0):
        return 3
    elif(x<120.0/365.0):
        return 1
    elif(x<125.0/365.0):
        return 2
    elif(x<170.0/365.0):
        return 1
    elif(x<228.0/365.0):
        return 2
    elif(x<302.0/365.0):
        return 1
    elif(x<=1.0):
        return 3
    else:
        print("time frac > 1. Returning full power")
        return 1

