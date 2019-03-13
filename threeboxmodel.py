import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#parameters
def rhs(s,T,*args):
    
    Fl,Fo,lambda_l,lambda_o,alpha_l,alpha_o,Cl,Co,Cdo,fl,gamma_o=args
    
    return [Fl[int(s)]/Cl-lambda_l/Cl*T[0]+alpha_o/(Cl*fl)*T[1]-alpha_l/(Cl*fl)*T[0],1/Co*(Fo[int(s)]-lambda_o*T[1]-alpha_o/(1-fl)*T[1]+alpha_l/(1-fl)*T[0]-gamma_o*(T[1]-T[2])),gamma_o/Co*(T[1]-T[2])]

from scipy.integrate import solve_ivp



Fl0=7.6
Fo0=6.5
#Fo0=0.
#Fl0=0.0

forcing_type="double_then_constant"
#forcing_type="1pctCO2"
#forcing_type="abrupt4xCO2"

if forcing_type=="double_then_constant":
#CO2 doubling, then held constant
    Fl=np.append(Fl0/2.*np.log(1.01)/np.log(2.)*np.arange(0,70),Fl0/2*np.ones_like(np.arange(70,141)))

    Fo=np.append(Fo0/2.*np.log(1.01)/np.log(2.)*np.arange(0,70),Fo0/2*np.ones_like(np.arange(70,141)))

elif forcing_type=="1pctCO2":
#1% until quadrupling
    Fl=7.6*np.log(1.01)/np.log(2.)*np.arange(0,140+1)
    Fo=6.5*np.log(1.01)/np.log(2.)*np.arange(0,140+1)
elif forcing_type=="abrupt4xCO2":
    Fl=Fl0*np.ones_like(np.arange(0,140+1))
    Fo=Fo0*np.ones_like(np.arange(0,140+1))

    

lambda_l=1.13

lambda_o=1.11


alpha_l=1.62
alpha_o=2.39
#alpha_l=alpha_o
Cl=1.
Co=12.
Cdo=244.
fl=0.33
gamma_o=2.04
params=[Fl,Fo,lambda_l,lambda_o,alpha_l,alpha_o,Cl,Co,Cdo,fl,gamma_o]


initial_value=[0,0,0]
interval=(0,140)
res=solve_ivp(fun=lambda s,T: rhs(s,T,*params),t_span=interval,y0=initial_value,t_eval=np.arange(0,140+1))
Tl,To,Tdo=res.y
Tdict={}
Tdict["Land"]=Tl
Tdict["Ocean"]=To
Tdict["Deep Ocean"]=Tdo
Fav=fl*Fl+(1-fl)*Fo
Tav=fl*Tl+(1-fl)*To
Q=gamma_o/Co*(To-Tdo)
lambda_av=fl*lambda_l+(1-fl)*lambda_o
lambda_inferred=(Fav-Q)[10:]/Tav[10:]


def colors(label):
    d={}
    d["Land"]=cm.Greens(.8)
    d["Ocean"]=cm.Blues(.4)
    d["Deep Ocean"]=cm.Blues(.9)
    d["abrupt4xCO2"]=cm.magma(.8)
    d["1pctCO2"]=cm.magma(.1)
    d["double_then_constant"]=cm.magma(.5)
    return d[label]

def plot_T(Tdict):
    for k in Tdict.keys():
        plt.plot(Tdict[k],label=k,color=colors(k))
        plt.legend(loc=0)


class SimpleModel():
    def __init__(self,forcing_type,num_years=140,lambda_l=1.13,lambda_o=1.11,alpha_l=1.62,alpha_o=2.39,Cl=1.,Co=12.,Cdo=244.,fl=0.33,gamma_o=2.04,Fl0=7.6,Fo0=6.5):
        if forcing_type=="double_then_constant":
            #CO2 doubling, then held constant
            Fl=np.append(Fl0/2.*np.log(1.01)/np.log(2.)*np.arange(0,num_years/2),Fl0/2*np.ones_like(np.arange(num_years/2,num_years+1)))

            Fo=np.append(Fo0/2.*np.log(1.01)/np.log(2.)*np.arange(0,num_years/2),Fo0/2*np.ones_like(np.arange(num_years/2,num_years+1)))


        elif forcing_type=="1pctCO2":
            #1% until quadrupling
            Fl=Fl0/2.*np.log(1.01)/np.log(2.)*np.arange(0,num_years+1)
            Fo=Fo0/2.*np.log(1.01)/np.log(2.)*np.arange(0,num_years+1)
        elif forcing_type=="abrupt4xCO2":
            Fl=Fl0*np.ones_like(np.arange(0,num_years+1))
            Fo=Fo0*np.ones_like(np.arange(0,num_years+1))

        params=[Fl,Fo,lambda_l,lambda_o,alpha_l,alpha_o,Cl,Co,Cdo,fl,gamma_o]
        
        initial_value=[0,0,0]
        interval=(0,num_years)
        res=solve_ivp(fun=lambda s,T: rhs(s,T,*params),t_span=interval,y0=initial_value,t_eval=np.arange(0,num_years+1))
        Tl,To,Tdo=res.y
        Tdict={}
        Tdict["Land"]=Tl
        Tdict["Ocean"]=To
        Tdict["Deep Ocean"]=Tdo
        self.Tdict=Tdict
        

        self.Fl=Fl
        self.Fo=Fo
        self.lambda_av=fl*lambda_l+(1-fl)*lambda_o
        
        self.Fav=fl*Fl+(1-fl)*Fo
        self.Tav=fl*Tl+(1-fl)*To
        self.Q=gamma_o/Co*(To-Tdo)
        
        self.lambda_inferred=(self.Fav-self.Q)[10:]/self.Tav[10:]

        self.num_years=num_years
    
    def plot_T(self,**kwargs):
        for k in self.Tdict.keys():
            plt.plot(self.Tdict[k],label=k,color=colors(k),**kwargs)
            plt.legend(loc=0)
    def plot_efficacy(self,**kwargs):
        plt.plot(np.arange(self.num_years+1)[10:],self.lambda_inferred/self.lambda_av,**kwargs)
        plt.ylabel(r'$\frac{\lambda_{inferred}}{\lambda}$')
        plt.xlabel("Years after start")
    def plot_efficacy_2(self):
        num=self.Tav/(self.Fav-self.Q)
    def plot_ECS(self,**kwargs):
        F2CO2=3.4482131495571124
        plt.plot(np.arange(self.num_years+1)[10:],F2CO2/self.lambda_inferred,**kwargs)
        plt.axhline(F2CO2/self.lambda_av)
       
        
    
OnePct=SimpleModel("1pctCO2")
DoublingThenConstant=SimpleModel("double_then_constant")
Abrupt=SimpleModel('abrupt4xCO2')


def solve_for_ECS(lambda_l=1.13,lambda_o=1.11,alpha_l=1.62,alpha_o=2.39,Cl=1.,Co=12.,Cdo=244.,fl=0.33,gamma_o=2.04,Fl0=7.6,Fo0=6.5):
    OnePct=SimpleModel("1pctCO2",lambda_l=lambda_l,lambda_o=lambda_o,alpha_l=alpha_l,alpha_o=alpha_o,Cl=Cl,Co=Co,Cdo=Cdo,fl=fl,gamma_o=gamma_o,Fl0=Fl0,Fo0=Fo0)
    DoublingThenConstant=SimpleModel("double_then_constant",lambda_l=lambda_l,lambda_o=lambda_o,alpha_l=alpha_l,alpha_o=alpha_o,Cl=Cl,Co=Co,Cdo=Cdo,fl=fl,gamma_o=gamma_o,Fl0=Fl0,Fo0=Fo0)
    Abrupt=SimpleModel('abrupt4xCO2',lambda_l=lambda_l,lambda_o=lambda_o,alpha_l=alpha_l,alpha_o=alpha_o,Cl=Cl,Co=Co,Cdo=Cdo,fl=fl,gamma_o=gamma_o,Fl0=Fl0,Fo0=Fo0)
    OnePct.plot_ECS(color=colors("1pctCO2"),label="1pctCO2")
    Abrupt.plot_ECS(color=colors('abrupt4xCO2'),label="abrupt4xCO2")
    DoublingThenConstant.plot_ECS(color=colors("double_then_constant"),label="double_then_constant")
    plt.legend()
    
    
                        
#(lambda_l=1.13,lambda_o=1.11,alpha_l=1.62,alpha_o=2.39,Cl=1.,Co=12.,Cdo=244.,fl=0.33,gamma_o=2.04,Fl0=7.6,Fo0=6.5):
