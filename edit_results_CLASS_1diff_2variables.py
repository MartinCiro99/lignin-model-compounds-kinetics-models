import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import os
from sklearn.metrics import r2_score
import plotly.graph_objects as go
from plotly.offline import plot as plott
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, curve_fit
import GPyOpt
from GPyOpt.methods import BayesianOptimization
import parr_editor


def find_reaction_file (reaction_no, print_search=False):
    for i in os.listdir(path):
        if print_search:
            print( f"i: {i}    i[: len(reaction_no)]: {i[: len(reaction_no)]}   {i[len(reaction_no)+1]}" )
        if i[: len(reaction_no)] == reaction_no and i[len(reaction_no)]== " ":
            break
    return i

def read_reaction():
    file=find_reaction_file(reaction_no, print_search=False) #string of filename
    file_path=os.path.join(path,file)
    df= pd.read_excel(file_path, index_col=0, skiprows=3)
    df.replace("n.a.", 0.0, inplace=True)
    df.rename(columns={df.columns[0]:"time"}, inplace = True)
    try:
        df["time"].iloc[1:-1]=df["time"].iloc[1:-1].apply(lambda x: x.split()[0]).astype(float)
    except:
        print("couldnt split times")
    return df

def calibration_equation(x,slope,intercept):
    y=slope*x+intercept
    return y
    
def read_calibration(reaction_no, component, plot=True):
    file=f"{cal_curve} {component} calibration.xlsx" #string of filename
    file_path=os.path.join(path,file)
    #df_cal= pd.read_excel(file_path,skiprows=5)
    df_cal= pd.read_excel(file_path,skiprows=3, names=["exp_no", "c","t_ret","area"])
    try:
        df_cal["c"]=df_cal["c"].apply(lambda x: x.split()[1]).astype(float)
    except:
        df_cal["c"].astype(float)
    
    slope,intercept = np.polyfit(x=df_cal.c , y=df_cal.area, deg=1)
    
    if plot:
        plt.clf()
        plt.rcParams["figure.figsize"]=(6,4)
        plt.rcParams["figure.dpi"]=300
        plt.rcParams["axes.labelsize"]=12
        #plt.rcParams["font.size"]=10
        x = np.linspace(0,df_cal.c.max(),100)
        y = calibration_equation(x=x, slope=slope, intercept=intercept)
        y_pred = calibration_equation(x=df_cal.c, slope=slope, intercept=intercept)
        r2= r2_score(y_true= df_cal.area, y_pred=y_pred)
        plt.scatter(df_cal["c"],df_cal["area"])
        plt.plot(x,y, "--", alpha=0.8)
        plt.text(x=df_cal.c.max()*0.0 , y=df_cal.area.max()*0.8, fontsize=12, s=f"$Area={round(slope,2)} \cdot C + {round(intercept,2)},$  $R^{{2}}={round(r2,4)}$")
        plt.xlabel("C [mg/ml]",fontstyle="italic")
        plt.ylabel("Area",fontstyle="italic")
        plt.title(f"{component} Calibration Curve")
        plt.savefig(f"{component}_calibration",dpi=300)
        plt.show();
    
    return df_cal,slope, intercept

def A_c (A): 
    return (A-intercept)/slope *dilution*cal_inj_vol/inj_vol/M[component]*1000

def create_reactants_df():
    parr_name= parr_editor.find_parr_name(reaction_no)
    df_parr = parr_editor.read_parr_file(parr_name=parr_name)
    df_parr_sampling = parr_editor.sampling_slice(df_parr=df_parr, T=T)
    starting_time_act, sampling_P = parr_editor.reaction_start_time_and_pressure(start_time, df_parr_sampling)
    print ("starting_time_act:", starting_time_act, " sampling_P: ", sampling_P, "time_dtype:", type(starting_time_act))
    
    reactants= df[df["time"].apply(lambda x: type(x) == float)].iloc[:,[0,2]]
    reactants[f"{component} [mmol/l]"] = reactants.iloc[:,1].apply(lambda A: A_c(A) )
    if reaction_no=="MC013":
        reactants["PO2"]=np.array([11.78,11.52,11.22,10.93,10.63,10.28,10.01,9.73,9.93,9.10,8.88])
    else:
        reactants["PO2"]= reactants["time"].apply(lambda i: parr_editor.reaction_start_time_and_pressure(
            starting_time_act+pd.Timedelta(minutes=i), df_parr_sampling)[1])
    reactants["O2 [mmol/l]"] = reactants["PO2"].apply(lambda i: C_pT(p=i,T=T))
    if reaction_no in ["MC025", "MC026", "MC027" ]:
        file=f"{reaction_no} {component} pH.xlsx"
        file_path=os.path.join(path,file)
        pHs=pd.read_excel(file_path)
        reactants.reset_index(inplace=True,drop=True)
        reactants["pH"]=pHs["pH"]
        reactants["OH"] = (10**-(14-reactants["pH"]))*1000 #mmol/l
    reactants.drop(columns=[reactants.columns.values[1]], inplace=True)
    reactants.reset_index(inplace=True,drop=True)
    if reaction_no in ["MC012"]:
        reactants.drop([7,9],axis=0,inplace=True)
    elif reaction_no =="MC017":
        reactants.drop([6],axis=0,inplace=True)
    return reactants

def plot_converted(save_fig=False):
    convertion=(1-(reactants.iloc[-1,1]/A_c( df.iloc[0,2])))*100 #.iloc ima vrstica,stolpec
    plt.clf()
    plt.plot(reactants["time"], reactants.iloc[:,1], "o")
    plt.axvline(x=0, color='grey',alpha=0.7,linestyle="-",lw=0.5)
    plt.plot(-10,A_c(df.iloc[0,2]),"ro")
    #plt.rcParams['font.family'] = 'sans-serif'
   #plt.rcParams['font.sans-serif'] = ['sans-serif']
    
    plt.ylabel(reactants.columns[1], fontstyle="italic")
    plt.xlabel("t [min]", fontstyle="italic")
    plt.ylim(top=df_cal["c"].max()/M[component]*1000*dilution*1.05, bottom=0)
    #plt.ylim(top=reactants.iloc[:,2].max()*1.2, bottom=0)
    plt.text(x=100 , y=reactants.iloc[:,2].max()*0.8,  s=f"convertion= {round(convertion,1)}%  $pH_{{final}}={pH_final}$")
    plt.title(f"{component}  T={T} °C  pH={pH}")
    plt.show()
    if save_fig:
        plt.savefig(f'{component}_{T}_{pH}_degradation.png', dpi=300)
    return

def C_pT(p,T): #p [bar], T [°C]
    T=T+273.15 #to K
    p=p*0.986923 #to atm
    C=p*np.exp((0.046*T**2 +203.357*T*np.log(T/298) -(299.378+0.092*T)* (T-298)- 20591)/ (8.3144*T))*1000
    return C #mmol/l

def plot_C_pT():
    p_values = np.linspace(1, 12, 100)
    T_values = np.linspace(20, 130, 100)
    p_grid, T_grid = np.meshgrid(p_values, T_values)
    C_grid = C_pT(p_grid, T_grid)
    fig = go.Figure(data=[go.Surface(
    z=C_grid, x=p_grid, y=T_grid,
    contours={
        "x": {"show": True, "start": 1, "end": 10, "size": 1, "color": "white"},
        "y": {"show": True, "start": 20, "end": 130, "size": 10, "color": "white"},
        "z": {"show": True, "start": np.min(C_grid), "end": np.max(C_grid), "size": 500, "color": "white"}})])
    fig.update_layout(
        title='Topnost kisika [mmol/ml]',
        scene=dict(
            xaxis_title='Pressure (bar)',
            yaxis_title='Temperature (°C)',
            zaxis_title='C (mmol/l)'))
    plott(fig, auto_open=True)
        
def plot_pH():
    plt.clf()
    plt.plot(reactants.time, reactants.pH,"o",reactants.time, reactants.OH,"o",reactants.time, reactants.iloc[:,1],"o")
    plt.legend(["pH","OH-","resorcinol"])
    return


def pH_fit(plot=True):
    if reaction_no in ["MC025","MC026","MC027"]: #resorcinol
        X_data=reactants.time.to_numpy()
        Y_data=reactants.pH.to_numpy()
    elif reaction_no in ["MC012"]:
        X_data=np.array([0,15,360])
        Y_data=np.array([pH_0,5.05,pH_final])
    elif reaction_no in ["MC015", "MC024"]:
        X_data=np.array([0,360])
        Y_data=np.array([5.05,pH_final])
    else:
        X_data=np.array([0,360])
        Y_data=np.array([pH_0,pH_final])
    #Y_data = 10**-(14-Y_data)
    initial_guesses = [10,0]
    params, covariance = curve_fit(pH_fit_model, X_data, Y_data, p0=initial_guesses)
    a_opt, b_opt = params[0], params[1]
    if plot:
        X_fit=np.linspace(0,360,360)
        Y_fit = pH_fit_model(X_fit, a_opt, b_opt)
        plt.scatter(X_data, Y_data, label='Experimental Data')
        plt.plot(X_fit, Y_fit,"--", label='Fitted Curve', color='red')
        plt.xlabel('t [min]',fontstyle="italic")
        plt.ylabel('pH', fontstyle="italic")
        txt_y=Y_data[0]-0.3*(Y_data[0]-pH_final)
        plt.text(x=110 , y=txt_y,fontsize=12,  s=f"$pH = {round(a_opt,4)}\cdot exp(-{round(b_opt,6)}\cdot t$)", fontstyle="italic")
        plt.title(f'{component}  T={T} °C')
        plt.show()
    print(f"Optimized Parameters:\n a = {a_opt}  b={b_opt}")
    return a_opt, b_opt
    
def pH_fit_model(X,a,b):
    #conv=(1-(reactants.iloc[-1,1]/df.iloc[1,1]))*100
    #OH_0,OH_f=10**-(14-pH_0),10**-(14-pH_final) 
    #dOH=abs(OH_0-OH_f)
    return a*np.exp(-X*b) #OH_0*np.exp(-X*a*dOH/(1-conv))

def activation_energy(k11,k12,k13,  k21,k22,k23,   T1=80,T2=100,T3=120,plot=True):
    R=8.314 #J/mol*K
    T1,T2,T3= 273.15+T1, 273.15+T2, 273.15+T3
    k11,k12,k13=np.log(k11),np.log(k12),np.log(k13)
    k21,k22,k23=np.log(k21),np.log(k22),np.log(k23)
    slp1,int1=np.polyfit(x=np.array([1000/T1,1000/T2,1000/T3]) ,y=np.array([k11,k12,k13]),deg=1)
    slp2,int2=np.polyfit(x=np.array([1000/T1,1000/T2,1000/T3]) ,y=np.array([k21,k22,k23]),deg=1)
    Ea1, Ea2 = (-slp1*R), (-slp2*R) #kJ/mol
    if plot:
        x=np.linspace(1000/T3, 1000/T1,100)
        #err1= r2_score(y_true=np.array([k11,k12,k13]) , y_pred=int1+np.array([k11,k12,k13])*slp1)
        #err2= r2_score(y_true=np.array([k21,k22,k23]) , y_pred=int2+np.array([k21,k22,k23])*slp2)
        
        plt.plot(np.array([1000/T1,1000/T2,1000/T3]) ,np.array([k11,k12,k13]), "o")
        plt.plot(x, int1+slp1*x, "--")
        plt.ylabel("$ ln(k_{1}) $", fontstyle="italic")
        plt.xlabel("1000/T [1/K]", fontstyle="italic")
        plt.text(x=1000/T2, y=k13, s= f"$Ea_{{1}}$ = {round(Ea1,4)} kJ/mol",fontsize=12)
        plt.savefig(f"{component} Ea1", dpi=300)
        plt.show()
        plt.clf()
        plt.plot(np.array([1000/T1,1000/T2,1000/T3]) ,np.array([k21,k22,k23]), "o")
        plt.plot(x, int2+slp2*x, "--")
        plt.ylabel("$ ln(k_{2}) $", fontstyle="italic")
        plt.xlabel("1000/T [1/K]", fontstyle="italic")
        plt.text(x=1000/T2, y=k23, s= f"$Ea_{{2}}$ = {round(Ea2,4)} kJ/mol",fontsize=12)
        plt.savefig(f"{component} Ea2", dpi=300)
        
    return Ea1, Ea2

def plot_LSE(opt_instance,mesh=False):
    fig = plt.figure(figsize=(14,12),dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(k_grid, m_grid, lse_grid, cmap='viridis')
    ax.scatter(ks,ms,LSEs,cmap=("viridis"))
    ax.set_xlabel('k')
    ax.set_ylabel('m')
    ax.set_zlabel('LSE')
    ax.set_zlim(0,35)
    ax.view_init(15, -180)
    plt.show()
    if mesh:
        k_values = np.linspace(0, 0.1, 10)
        m_values = np.linspace(0, 0.02, 10)
        k_grid, m_grid = np.meshgrid(k_values, m_values)
        lse_grid = np.zeros((len(k_values),len(m_values)))
    
        for i in range(len(k_values)):
            for j in range(len(m_values)):
                #print(i,j)
                lse_grid[i, j] = opt_instance.model_scipy_minimize(np.array([k_values[i],m_values[j],1]))  # Assuming n = 0
    
        fig = plt.figure(figsize=(14,12),dpi=300)
        ax = fig.add_subplot(111, projection='3d')
        #ax.plot_surface(k_grid, m_grid, lse_grid, cmap='viridis')
        ax.scatter(k_grid,m_grid,lse_grid,cmap=("viridis"))
        ax.set_xlabel('k')
        ax.set_ylabel('m')
        ax.set_zlabel('LSE')
        ax.set_zlim(0,35)
        ax.view_init(15, -180)
        plt.show()
        return
    
        
class opt:
    def __init__(self,k,m=1,n=1):
        self.k=k
        self.m=m
        self.n=n
    def sistem_ode_resorcinol(self,c,t):
        A=c[0]
        B=c[1]
        #C=c[2]
        if np.argmax(t < reactants.iloc[:, 0].values) == self.sampled:
            self.sampled+=1
            #self.C+=-0.3
            self.C=self.pressures[self.sampled-1]
            #print("sampled, C:",self.C)
        if self.optimize ==1:
            dAdt= -self.k* A * (1000*Kw*10**B) * self.C
            dBdt=-a*b*np.exp(-b*t)
        elif self.optimize ==2:
            H=1000*10**(-B)
            #dAdt= -self.k* A * (1000*Kw*10**B) * self.C
            dAdt = - (Ka1* (self.k*1000*10**(-B)+self.m*Ka2)) / (1000*10**(-2*B)+Ka1*(1000*10**-B)+Ka1*Ka2) *A *self.C
            #dAdt = -self.k* A**self.m *  (1000*Kw*10**B)**self.n *self.C
            #dAdt=-self.k* A**self.m * (1000*Kw*10**B)**self.n  * self.C
            #dAdt= -((self.k/(1+Ka1*H)) + (self.m*Ka1/H))  *A**self.n*self.C ##vanillin article
            #dAdt=-self.k*(1-self.k/self.m)*(1-1/(1+Ka1/H))*A*self.C
            #dAdt=-self.k*(Kw/H)**self.m*self.C**self.n
            #dAdt=-(self.k*H/(Ka1+H)+self.m*Ka1/(Ka1+H))*A**1 *self.C**1
            #dAdt = -self.k*H/(Ka1+H) *A*self.C**3 -self.m*Ka1/(H+Ka1)*A*self.C**1
            #dAdt = -self.k* ((self.m*H)/(1+self.m*H))**2 * A**2
            #dAdt=-self.k*A**self.m*self.C**(3/2)
            dBdt=-a*b*np.exp(-b*t)
        elif self.optimize == 3:
            H=1000*10**(-B)
            #dAdt= -(self.k +self.m*Ka1/(10**-B) +self.n*Ka1*Ka2/(10**-2*B))/(1+ Ka1/(10**-B) + Ka1*Ka2/(10**-2*B)) * A * self.C
            dAdt=-self.k* A**self.m  *(Kw/H)**self.n *self.C**1
            #dAdt=-self.k* A**self.m  * self.C**self.n
            #dAdt=-self.k*self.m*((1+self.n*H)/(1+self.n/self.m*H))*A*self.C
            dBdt=-a*b*np.exp(-b*t)
        #print("A:",A,"  B:",B,"  C:",self.C,"  t:",t)

        dy1= dAdt
        dy2= dBdt
        return ([dy1,dy2])
    
    def resevanje_ode_resorcinol(self,k_opt=False):
        #print("2.res_ODE",self.k)
        self.steps=361
        self.C=reactants.iloc[0,3]
        self.pressures=reactants.iloc[:,3].to_numpy()
        self.sampled=1
        t = np.linspace(0.0,360.0,self.steps)
        if reaction_no in ["MC025","MC026","MC027"]:
            zacetni = reactants.iloc[0,[1,4]] #resorcinol, ph
        else:
            zacetni = [reactants.iloc[0,1],pH_0]
        #print("beginning odeint",self.sampled)
        resitev_modela = odeint(self.sistem_ode_resorcinol,zacetni,t=t)
        A,B = resitev_modela[:,0], resitev_modela[:,1]
        if k_opt:
            return A,B
        indices=reactants.time.values*(self.steps-1)/360. #prej bilo *self.steps/360. -1
        #print("pred")
        #print(A, len(A),type(A),A.shape,indices.astype(int).tolist())
        Ao = np.take(A,indices.astype(int).tolist())
        #print("po")
        #print("odeint opravljen")
        return Ao
    
    def scipy_fit(self,optimize=2, plot=True, scipy_curve_fit=False):
        self.optimize=optimize
        self.A_lab=reactants.iloc[:,1]
        initial_guesses=[0.3,2,0.5]
        bounds=[(0,3.51*10**(2)),(2,2),(0.5,0.5)] #[(1e-7,2),(0,1),(1,1)] for resorcinol
        X_data=reactants.time
        Y_data=reactants.iloc[:,1]
        result = minimize(self.model_scipy_minimize, initial_guesses, method="Nelder-Mead",bounds=bounds, options={'maxiter': 10000, 'gtol': 1e-8})
        opt_params = result.x #coordinates - k,m
        print("Optimized parameters:", opt_params)
        if plot:
            A,B=self.resevanje_ode_resorcinol(k_opt=True)
            t=np.linspace(0,360,self.steps)
            plt.clf();
            plt.axvline(x=0, color='grey',alpha=0.7,linestyle="-",lw=0.5)
            plt.plot(-10,A_c(df.iloc[0,2]),"ro")
            plt.rcParams['font.family'] = 'sans-serif'
            A_lab=reactants.iloc[:,1]
            plt.plot(reactants.time, A_lab,"o")
            plt.plot(t,A,"--")
            plt.xlabel("t [min]", fontstyle="italic")
            plt.ylabel("C [mmol/l]", fontstyle="italic")
            plt.title(f"{component}  T={T} °C  pH={pH}")
            #plt.ylim(0,22)
            #plt.text(x=75,y=25,fontsize=12 ,s=f"$k_{{1}}$={round(opt_params[0],5)}  $k_{{2}}$={round(opt_params[1],9)}") #   n={round(opt_params[2],6)}") :for resorcinol: f"$k_{{1}}$={round(opt_params[0],6)}  $k_{{2}}$={round(opt_params[1],4)}"    f"k={round(opt_params[0],7)} $l^{{2.5}}mmol^{{-2.5}} $"

            plt.text(x=75,y=17.5,fontsize=12 ,s=f"$k_{{1}}$={round(opt_params[0],6)}  $k_{{2}}$={round(opt_params[1],8)}") #   n={round(opt_params[2],6)}") :for resorcinol: f"$k_{{1}}$={round(opt_params[0],6)}  $k_{{2}}$={round(opt_params[1],4)}"    f"k={round(opt_params[0],7)} $l^{{2.5}}mmol^{{-2.5}} $"
            #plt.text(x=95,y=A_lab.max()*0.95,fontsize=12 ,s=f"$k_{{1}}$={round(opt_params[0],3)}  m={round(opt_params[1],3)}  n={round(opt_params[2],3)}")
            #"$k_{{1}}$={round(opt_params[0],9)}  $k_{{2}}$={round(opt_params[1],4)}"
            #f"$k_{{1}}$={round(opt_params[0],7)}  $k_{{2}}$={round(opt_params[1],3)}"
            #f"$k_{{1}}$={round(opt_params[0],8)}  m={round(opt_params[1],4)}  n={round(opt_params[2],4)}"
            plt.savefig(f"{component}_model_fit_{T}_{pH_0}.png",dpi=300)
            plt.show()
            LSE=self.model_scipy_minimize(np.array([self.k,self.m,self.n]))
            print(LSE)
        return result
        if scipy_curve_fit:
            params, covariance = curve_fit(self.model_curve_fit, X_data, Y_data, p0=initial_guesses,bounds=([0,0], [1,5]),maxfev=5000)
            k_opt, m_opt= params[0], params[1]
            print(k_opt, m_opt)
            if plot:
                A,B=self.resevanje_ode_resorcinol(k_opt=True)
                t=np.linspace(0,360,self.steps)
                plt.clf();
                A_lab=reactants.iloc[:,1]
                plt.plot(reactants.time, A_lab,"o")
                plt.plot(t,A,"--")
                plt.text(x=120,y=A_lab.max()*0.7,s=f"k={k_opt}  m={m_opt}")
                plt.show()
            return params, covariance
        
    def plot_model(self,x):
        self.k= x[0]
        self.m=x[1]
        self.n=x[2]
        LSE=self.model_scipy_minimize(x)
        print(LSE)
        A,B=self.resevanje_ode_resorcinol(k_opt=True)
        t=np.linspace(0,360,self.steps)
        plt.clf();
        plt.axvline(x=0, color='grey',alpha=0.7,linestyle="-",lw=0.5)
        plt.plot(-10,A_c(df.iloc[0,2]),"ro")
        A_lab=reactants.iloc[:,1]
        plt.plot(reactants.time, A_lab,"o")
        plt.plot(t,A,"--")
        plt.xlabel("t [min]", fontstyle="italic")
        plt.ylabel("C [mmol/l]", fontstyle="italic")
        plt.title(f"{component}  T={T} °C  pH={pH}")
        plt.text(x=100,y=A_lab.max()*0.95,fontsize=12 ,s=f"$k_{{1}}$={round(self.k,8)}  $k_{{2}}$={round(self.m,6)}") #   n={round(opt_params[2],6)}") :for resorcinol: f"$k_{{1}}$={round(opt_params[0],6)}  $k_{{2}}$={round(opt_params[1],4)}"    f"k={round(opt_params[0],7)} $l^{{2.5}}mmol^{{-2.5}} $"
        plt.show()
        return
        
        
    def model_scipy_minimize(self, x):
        self.k= x[0]
        #print(self.k)
        self.m=x[1]
        self.n=x[2]
        Ao=self.resevanje_ode_resorcinol()
        LSE= 1/len(Ao)*np.sum( (Ao - self.A_lab)**2) #least square erro
        #print(Ao, "\n", self.A_lab.to_numpy())
       # print("odeint done",self.k,self.m, LSE)
        ks.append(self.k), ms.append(self.m), LSEs.append(LSE)
        return LSE
        
    def model_curve_fit(self, A, k, m):
        self.k=k
        self.m=m
        Ao=self.resevanje_ode_resorcinol()
        LSE= 1/len(Ao)*np.sum( (Ao - self.A_lab)**2)
        print("LSE", LSE)
        return LSE #vrednost, katere minimum hočeš, tj. med A in A_lab
    
    ##for gaussian
    def k_LSE(self,x):
        self.x=x.reshape((self.optimize,))
        if self.optimize==1:
            try:
                self.k=10**(-x[0,0])
            except:
                self.k= 10**(-x)
        elif self.optimize==2:
            self.k= 10**(-self.x[0])
            self.m=self.x[1]
        print("k", self.k)
        Ao=self.resevanje_ode_resorcinol()
        A_lab=reactants.iloc[:,1]
        LSE= 1/len(Ao)*np.sum( (Ao-A_lab)**2)
        return  LSE
    
    def find_my_minimum(self, optimize=1):
        self.optimize=optimize
        if self.optimize ==1:
            domain=[{'name' : "var_1", 'type': 'continuous', 'domain': (1E-8,8)}]
        elif self.optimize ==2:
            domain=[{'name' : "var_1", 'type': 'continuous', 'domain': (1E-8,8)},
                    {'name' : "var_2", 'type': 'continuous', 'domain': (0,2)}]
        myBopt_1d = BayesianOptimization(f=lambda x: self.k_LSE(x), domain=domain,
                                         model_type="GP",acquisition_type="EI")
        myBopt_1d.run_optimization(max_iter=50, verbosity=True,)
        GPY_results=myBopt_1d
        print("sampled:", self.sampled, "C:",self.C)
        return GPY_results
    
    def plot_GPY_model(self):
        plt.clf()
        self.k=10**-GPY_results.x_opt[0]
        A,B=self.resevanje_ode_resorcinol(k_opt=True)
        plt.plot(reactants.time,reactants.iloc[:,1],"o", np.linspace(0.0,360.0,self.steps),A,"--")
        if self.optimize==2:
            self.m=GPY_results.x_opt[1]
        plt.text(x=0,y=reactants.iloc[1,1]*0.8, s=f"k={self.k}  m={self.m}")    
        plt.show()
        return A,B
    def plot_LSE(self):
        plt.clf()
        if self.optimize==1:
            plt.plot(ks,diff,"o")
            plt.show()
        elif self.optimize==2:
            k,m=ks[:,0], ks[:,1]
            LSE=diff
            fig = plt.figure(figsize=[14,12],dpi=300)
            ax = fig.add_subplot(111, projection='3d')
            p=ax.scatter(k, m, LSE, c=LSE, marker='o',cmap=plt.cm.viridis)
            ax.set_zlim(0,15)
            # Set labels
            #ax.title("Least square error for k1 and k2")
            ax.set_xlabel('log(k1)')
            ax.set_ylabel('$ n $')
            ax.set_zlabel('LSE')
            ax.view_init(elev=15, azim=90)
            
            #norm = plt.Normalize(LSE.min(), LSE.max())
            #colors = cm.viridis(norm(LSE))
            #mappable = cm.ScalarMappable(norm=norm, cmap='viridis')
            #mappable.set_array(LSE)
            fig.colorbar(p, ax=ax, shrink=0.6, aspect=10)
            plt.show()


path="C:\\Users\\Admin\\OneDrive\\MAG1\\.raziskovalna\\.results\\reactions_results" #where to look for data
M={"Vanillin":      152.15, "Catechol":      110.1, "Resorcinol":    110.1, "Acetovanillone":166.174}
Ka1,Ka2,Kw=[4.8E-10, 7.9E-12, 1.0E-14] #resorcinol only
#Ka1,Ka2,Kw=[3.24E-10, 9.33E-14, 1.0E-14] #catechol
#Ka1,Kw = [4.0179E-8, 1.0E-14] #vanillin
#Ka1,Kw = [10**-8.17, 1.0E-14] #acetovanillone
k_init=0.0014
#-------
reaction_no="MC027"
#--------

Reactions_table = pd.read_excel(os.path.join("C:\\Users\\Admin\\OneDrive\\MAG1\\.raziskovalna\\.results","Reactions_table.xlsx"))
Reactions_table["start_time"] =  pd.to_datetime(Reactions_table["start_time"], format='%H:%M:%S')
component, T, pH, pH_0, pH_final, dilution, inj_vol, cal_curve, P_cal_curve, cal_inj_vol, start_time = Reactions_table[Reactions_table["reaction_no"] == reaction_no].iloc[:,1:].values[0,:]
df_cal, slope, intercept = read_calibration(reaction_no=cal_curve,component=component, plot=False)
df = read_reaction()
reactants= create_reactants_df()
plot_converted(save_fig=False)
"""
#plot_C_pT()
#plot_pH()
a,b=pH_fit()

convertion=1-(reactants.iloc[:,1]/A_c(df.iloc[0,2]))
opt_instance = opt(k_init)
ks,ms,LSEs=[],[],[]
result=opt_instance.scipy_fit(optimize=2)

k13,k23=result.x[0],result.x[1]"""
#opt_instance.model_scipy_minimize(np.array([0.0077,0,1]))
#opt_instance.plot_model(np.array([10**(-3.45),2*10**(-2),1]))
#plt.plot(convertion, reactants.OH,"o")

#gaussian optimisation
""""
GPY_results=opt_instance.find_my_minimum(optimize=2)
ks=GPY_results.X
diff=GPY_results.Y
x_opt=GPY_results.x_opt
A,B=opt_instance.plot_GPY_model()
opt_instance.plot_LSE()
GPY_results.x_opt
"""
