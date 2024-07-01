import pandas as pd
import numpy as np
import re
import os 
#first= pd.read_csv("MC013_Acetovanillone_120-8-10__first.csv", sep="(.*?) (.*?),(\d{1,2},\d{1,2}),(\d{1,2},\d{1,2}),(\d{1,2},\d{1,2}),(-?\d{1,2},\d{1,2})",engine="python",skiprows=0)
#second= pd.read_csv("MC013_Acetovanillone_120-8-10__second.csv", sep=",")
#second.reset_index(inplace=True)
#first.reset_index(inplace=True)

parr_name = "MC027_Resorcinol_65-12-10.csv"
parr_path = "C:\\Users\\Admin\\OneDrive\\MAG1\\.raziskovalna\\.results\\Parr logs"
broken_names=["MC013_Acetovanillone_120-8-10__first.csv", "MC013_Acetovanillone_120-8-10__second.csv"]


def read_parr_file(parr_name):
    if parr_name==["MC013_Acetovanillone_120-8-10__first.csv", "MC013_Acetovanillone_120-8-10__second.csv"]:
        df_concated=broken_parr(["MC013_Acetovanillone_120-8-10__first.csv", "MC013_Acetovanillone_120-8-10__second.csv"])
        return df_concated
    f=open(os.path.join(parr_path, parr_name), "r")
    f.readline()
    f=f.read()
    matches = re.findall(r"(.*?) (.*?),(\d{1,3},\d{1,3}),(\d{1,3},\d{1,2}),(\d{1,2},\d{1,2}),(-?\d{1,2},\d{1,2})", f)
    columns=["date","time", "T", "T_set", "P", "dp"]
    df=pd.DataFrame(data=matches, columns=columns)
    df.drop(labels=["T_set"], axis=1,inplace=True)
    df.replace(",",".", inplace=True)
    if len(matches) == 0:
        matches = [tuple(re.findall(r'[^,\s]+', line)) for line in f.split('\n')]
        try:
            if "MC027" in parr_name:
                columns=["date","time", "T","stirring", "T_set", "P", "dp"]
            elif "MC026" in parr_name:
                columns=["date","time", "T","stirring", "T_set", "P", "dp"]
            else:
                columns=["date","time", "T", "T_set", "P", "dp","stirring"]
        
            df=pd.DataFrame(data=matches, columns=columns)
            df.drop(labels=["T_set","stirring"], axis=1,inplace=True)
        except:
            columns=["date","time", "T", "T_set", "P", "dp"]
            df=pd.DataFrame(data=matches, columns=columns)
            df.drop(labels=["T_set"], axis=1,inplace=True)
        #matches = re.findall(r"[^\s,]+", f)
    df["T"] = df["T"].str.replace(",",".").astype(float)
    df["P"] = df["P"].str.replace(",",".").astype(float)
    df["dp"] = df["dp"].str.replace(",",".").astype(float)
    df['time'] = pd.to_datetime(df['time'], format='%H:%M:%S')
    #df['time'] = df['time'].dt.strftime('%H:%M:%S')
    return df

def broken_parr(broken_names):
    for i in range(len(broken_names)):
        f=open(os.path.join(parr_path ,broken_names[i]) , "r")
        f.readline()
        f=f.read()
        matches = re.findall(r"(.*?) (.*?),(\d{3}),(\d{1,2},\d{3,4}),(\d{2,3},\d{3,4}),(-?\d{1,2},\d{3,4}),(\d{1},\d{3,4})", f)
        columns=["date","time", "num", "P", "T", "dp", "num2"]
        if i==0:
            df0=pd.DataFrame(data=matches, columns=columns)
            df0.drop(["num","num2"],axis=1, inplace=True)
        if i==1:
            df1=pd.DataFrame(data=matches, columns=columns)
            df1.drop(["num","num2"],axis=1, inplace=True)
            df_concated=pd.concat([df0, df1], join="outer",axis=0)
            df_concated.reset_index(inplace=True, drop=True)
            df_concated=df_concated.reindex(columns=["date","time","T","P","dp"])
            df_concated["T"] = df_concated["T"].str.replace(",",".").astype(float)
            df_concated["P"] = df_concated["P"].str.replace(",",".").astype(float)
            df_concated["dp"] = df_concated["dp"].str.replace(",",".").astype(float)
            df_concated['time'] = pd.to_datetime(df_concated['time'], format='%H:%M:%S')
            #df_concated['time'] = df_concated['time'].dt.strftime('%H:%M:%S')
    return df_concated

def find_parr_name(reaction_no="MC013"):
    if reaction_no=="MC013":
        return ["MC013_Acetovanillone_120-8-10__first.csv", "MC013_Acetovanillone_120-8-10__second.csv"]
    for i in os.listdir(parr_path):
        #print( f"i: {i}    i[: len(reaction_no)]: {i[: len(reaction_no)]}   {i[len(reaction_no)+1]}" )
        if i[: len(reaction_no)] == reaction_no and i[len(reaction_no)]== "_":
            break
    return i
    
#df_parr=read_parr_file(parr_name)
#df333= broken_parr(broken_names)

# finding the sampling time: slicing with stable T condition and dp boundaries
def sampling_slice(df_parr, T):
    #df_parr_sampling=df_parr[(df_parr["dp"]>-0.5)& (df_parr["dp"]<0.0) & (df_parr["T"]>T-3)& (df_parr["T"]<T+2)]
    df_parr_sampling=df_parr[(df_parr["T"]>T-3)& (df_parr["T"]<T+2)]
    return df_parr_sampling

####find starting time
#start_time=pd.Timestamp("01-01-1900 9:35:00") ##that is read from the Reactions_table.xlsx file

def reaction_start_time_and_pressure(start_time, df_parr_sampling, dp02=False):
    taking_sample= df_parr_sampling[(df_parr_sampling["time"]<start_time+pd.Timedelta(minutes=3)) & (df_parr_sampling["time"]> start_time-pd.Timedelta(minutes=3))]
    dp_min=taking_sample.dp.min()
    start_time_act= taking_sample[taking_sample.dp==dp_min]["time"].iloc[0]
    sampling_P=df_parr_sampling[( df_parr_sampling["time"] <=(start_time_act - pd.Timedelta(seconds=60)) ) & 
                       ( df_parr_sampling["time"] >(start_time_act - pd.Timedelta(seconds=90)))]["P"].mean()
    return start_time_act, sampling_P
#df=read_parr_file(parr_name)
#df333= broken_parr(broken_names)
#df_parr_sampling = sampling_slice(df_parr=df333, T=120)
parr_name= find_parr_name(reaction_no="MC013")
df_parr = read_parr_file(parr_name=parr_name)
df_parr_sampling=sampling_slice(df_parr,T=120)


