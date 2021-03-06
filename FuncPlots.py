import matplotlib.pyplot as plt
import scipy # is needed so.. just do it
import numpy as np

# function used for plotting
class FigSetting:
    FigTitle=""
    Title=""
    xlabel=""
    ylabel=""
    label=[]
    figsize=(10,7)
    FontsizeTitle=18
    Fontsize_label=15
    SetGrid=1
    SetTick=1
    SetLegend=1
    
    def __init__(self,_FigTitle="TitleFigure",_Title="Title",_xlabel="x label",_ylabel="y label",_label=[]):
        self.FigTitle=_FigTitle
        self.Title=_Title
        self.xlabel=_xlabel
        self.ylabel=_ylabel
        self.label=_label
    
    def Flabel(self,i):
        if(i>=len(self.label)):
            return "Species "+str(i)
        else:
            return self.label[i]
        

def FigPlotTime(FS,Sollution):# where solution is the return of solve_ivp
    # do event handeling
    if(type(Sollution)!=scipy.integrate._ivp.ivp.OdeResult):
        print("Error: invalid input Sollution. Given type: '"+type(Sollution).__name__+"' While expected: 'scipy.integrate._ivp.ivp.OdeResult'")
        return 0
    if(not Sollution.success):
        print("Error: Given input is of a failed intergration. Pls provide a solution of an intergration that was succesfull")
        return 0
    if(type(FS)!=FigSetting):
        print("Error: invalid input FS. Given type:'"+ type(FS).__name__+"' While expected: 'FigSetting'")
        return 
    # create figure and axes
    fig=plt.figure(FS.FigTitle, figsize=FS.figsize)
    s=fig.add_subplot(111)
    
    # Get the correct settings
    s.set_title(FS.Title,fontsize = FS.FontsizeTitle)
    s.set_xlabel(FS.xlabel,fontsize=FS.Fontsize_label)
    s.set_ylabel(FS.ylabel,fontsize=FS.Fontsize_label)
    if(FS.SetTick): s.tick_params(axis='both', labelsize=15)
    if(FS.SetGrid): s.grid()
    
    # plot the solution on to the figure
    for i in range(Sollution.y.shape[0]):
        s.plot(Sollution.t,Sollution.y[i],label=FS.Flabel(i))
    if(FS.SetLegend): s.legend(prop={'size': 20})
    return fig

def FigPlotPhaseSpace(FS,Sollution,dXdt):
    # do event handeling
    if(type(Sollution)!=scipy.integrate._ivp.ivp.OdeResult):
        print("Error: invalid input Sollution. Given type: '"+type(Sollution).__name__+"' While expected: 'scipy.integrate._ivp.ivp.OdeResult'")
        return 0
    if(not Sollution.success):
        print("Error: Given input is of a failed intergration. Pls provide a solution of an intergration that was succesfull")
        return 0
    if(type(FS)!=FigSetting):
        print("Error: invalid input FS. Given type:'"+ type(FS).__name__+"' While expected: 'FigSetting'")
        return 
    # create figure and axes
    fig=plt.figure(FS.FigTitle, figsize=FS.figsize)
    s=fig.add_subplot(111)
    
    # Get the correct settings
    s.set_title(FS.Title,fontsize = FS.FontsizeTitle)
    s.set_xlabel(FS.xlabel,fontsize=FS.Fontsize_label)
    s.set_ylabel(FS.ylabel,fontsize=FS.Fontsize_label)
    if(FS.SetTick): s.tick_params(axis='both', labelsize=15)
    if(FS.SetGrid): s.grid()
    
    s.plot(Sollution.y[0],Sollution.y[1])
    Gridpoints=20
    x = np.linspace(s.get_xlim()[0], s.get_xlim()[1], Gridpoints)
    y = np.linspace(s.get_ylim()[0], s.get_ylim()[1], Gridpoints)
    
    X1 , Y1  = np.meshgrid(x, y) 
    DX1, DY1 = dXdt(0,[X1, Y1])
    
    M = (np.hypot(DX1, DY1))
    M[ M == 0] = 1. # no devision by zero
    DX1 /= M # normalise
    DY1 /= M
    
    Q = s.quiver(X1, Y1, DX1, DY1, M, pivot='mid', cmap=plt.cm.jet)
    
    
    return fig