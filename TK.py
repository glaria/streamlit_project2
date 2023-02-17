#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      aglaria
#
# Created:     18/03/2019
# Copyright:   (c) aglaria 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import tkinter as tk
from functools import partial
from Ztest import *
# import numpy as np

# INPUT
#uplift = 0.01
#N = 8000 # Total size
#significance = 0.1
minimum_aceptance_control_group = 0.00 # Final split must be significant for all values between minimum_aceptance_control_group & maximum_aceptance_control_group
#maximum_aceptance_control_group = 0.05

epsilon = 0.001
step = 0.005
initial_solution = 0.1 # Initial split (for CG from total)
k = [0] #to select different modes

def evaluate(value):
    for cg_acceptance in cg_acceptance_list:
        tg_acceptance = uplift + cg_acceptance
        tg_acceptors = tg_acceptance*(1-value)*N
        cg_acceptors = cg_acceptance*value*N
        tgsize = (1-value)*N
        cgsize = value*N
        if (z2p(zscore(tg_acceptance, cg_acceptance, tgsize, cgsize )) >= significance - epsilon): #and cg_acceptors < 1
            return False
    print(z2p(zscore(tg_acceptance, cg_acceptance, tgsize, cgsize )))
    if maximum_aceptance_control_group*value*N < 2: #minimum number of records
        return False
    else:
        return True

def descend(init_value):
    value = init_value - step
    aux_sol.append(value)
    if evaluate(value) and value > 0.02:
        print("Descending to find a better solution")
        descend(value)
        print(value)
    elif (value <= 0.02) or (not evaluate(value)):
        solution.append(aux_sol[len(aux_sol)-1])

def ascend(init_value):
    value = init_value + step
    if not evaluate(value) and value < 0.5:
        print("Ascending to find a solution")
        ascend(value)
    elif evaluate(value):
        solution.append(value)

def calculate():
    if evaluate(initial_solution):
        descend(initial_solution)
        #Show TG/CG size, TG/CG acceptors and p-value
        cg_size = int(solution[0]*N)
        tg_size = N-cg_size
        print("CG size",cg_size)
        print("TG size", tg_size)
        return (solution[0],tg_size,cg_size)
    else:
        ascend(initial_solution)
        if len(solution) > 0:
            #Show TG/CG size, TG/CG acceptors and p-value
            cg_size = round(solution[0]*N)
            tg_size = N-cg_size
            print("CG size",cg_size)
            print("TG size", tg_size)
            return (solution[0],tg_size,cg_size)

        else:
            print("No solution")
            return (0, N, 0)

def adapt_cg(input_cg_acceptance, total_counts):
    if float(total_counts)*input_cg_acceptance/2 >= 3:
        return input_cg_acceptance
    else:
        return 6/float(total_counts)



def new_mode(labelTitle, labelNum1, labelNum2, labelNum3, labelNum4, label_result,label_tg, label_cg, labelp, n1, n2, n3, n4):

    k.append(k[len(k)-1]+1)
    n1.set('')
    n2.set('')
    n3.set('')
    n4.set('')
    if k[len(k)-1]%2 == 0:
        root.title('TG/CG min split')
        labelTitle.config(text="Getting min CG size")
        labelNum1.config(text="Total counts:")
        labelNum2.config(text="Expected uplift:")
        labelNum3.config(text="Significance:")
        labelNum4.config(text="Max. acceptance in control:")
    else:
        root.title('Z-score calculator')
        root.geometry('700x200+100+200')
        labelTitle.config(text="Getting Z-score and p-value")
        labelNum1.config(text="Target Group counts:")
        labelNum2.config(text="Target Group acceptors:")
        labelNum3.config(text="Control Group counts:")
        labelNum4.config(text="Control Group acceptors:")
    #print k
    label_result.config(text="")
    label_tg.config(text="" )
    label_cg.config(text="")
    labelp.config(text="")


def call_result(label_result,label_tg, label_cg, labelp, n1, n2, n3, n4):
    try:
        if k[len(k)-1]%2 == 0:
            global solution
            solution = []
            global aux_sol
            aux_sol = []
            global N
            N = float(n1.get())
            global uplift
            if n2.get()[len(n2.get()) -1] == '%':
                uplift = float(n2.get()[:-1])/100
            else:
                uplift = float(n2.get())
            global significance
            if n3.get()[len(n3.get()) -1] == '%':
                significance = float(n3.get()[:-1])/100
            else:
                significance = float(n3.get())
            global maximum_aceptance_control_group
            #maximum_aceptance_control_group = adapt_cg(float(n4.get()),N)
            if n4.get()[len(n4.get()) -1] == '%':
                maximum_aceptance_control_group = float(n4.get()[:-1])/100
            else:
                maximum_aceptance_control_group = float(n4.get())
            global cg_acceptance_list
            cg_acceptance_list = np.arange(maximum_aceptance_control_group, minimum_aceptance_control_group, -0.001)

            result = calculate()
            print(result)
            if result[0] > 0:
                label_result.config(text="Result: %f" % float(result[0]))
                label_tg.config(text="TG size: %d" % float(result[1]))
                label_cg.config(text="CG size: %d" % float(result[2]))
                p_value = z2p(zscore(maximum_aceptance_control_group+uplift, maximum_aceptance_control_group,result[1], result[2] ))
                labelp.config(text= "P-Value: %f" %round(p_value, 4))
                print(p_value)

            else:
                label_result.config(text="No solution")
                label_tg.config(text="" )
                label_cg.config(text="")
                labelp.config(text="")
            return
        else:
            global N1
            N1 = int(n1.get())
            global N2
            N2 = int(n3.get())
            global acceptor_n1
            acceptor_n1 = int(n2.get())
            global acceptor_n2
            acceptor_n2 = int(n4.get())
            print([N1,N2,acceptor_n1, acceptor_n2, float(acceptor_n1)/float(N1)])
            p_value = z2p(zscore(float(acceptor_n1)/float(N1), float(acceptor_n2)/float(N2),float(N1), float(N2)))
            acceptance_tg = acceptor_n1/float(N1)
            acceptance_cg = acceptor_n2/float(N2)
            uplift = str(round((acceptance_tg - acceptance_cg)*100,3))
            labelp.config(text="P-Value: %s" %p_value )
            label_result.config(text="Uplift: %s %%" %uplift )
            label_tg.config(text="Acceptance TG: %s" %str(round(acceptance_tg,3)))
            label_cg.config(text="Acceptance CG: %s" %str(round(acceptance_cg,3)))
    except ValueError:
        label_result.config(text="Input Error")
        labelp.config(text="")
        label_tg.config(text="")
        label_cg.config(text="")

root = tk.Tk()
root.geometry('600x200+100+200')
root.title('TG/CG min split')
number1 = tk.StringVar()
number2 = tk.StringVar()
number3 = tk.StringVar()
number4 = tk.StringVar()

labelTitle = tk.Label(root)
labelTitle.grid(row=0, column=3)
labelTitle.config(text="Getting min CG size")
labelNum1 = tk.Label(root)
labelNum1.grid(row=1, column=0)
labelNum1.config(text="Total counts:")
labelNum2 = tk.Label(root)
labelNum2.grid(row=2, column=0)
labelNum2.config(text="Expected uplift:")
labelNum3 = tk.Label(root)
labelNum3.grid(row=1, column=4)
labelNum3.config(text="Significance:")
labelNum4 = tk.Label(root)
labelNum4.grid(row=2, column=4)
labelNum4.config(text="Max. acceptance in control:")


labelResult = tk.Label(root)
labelResult.grid(row=7, column=2)
labelTGsize = tk.Label(root)
labelTGsize.grid(row=8, column=2)
labelCGsize = tk.Label(root)
labelCGsize.grid(row=9, column=2)
labelpValue = tk.Label(root)
labelpValue.grid(row=8,column = 4)
entryNum1 = tk.Entry(root, textvariable=number1).grid(row=1, column=2)
entryNum2 = tk.Entry(root, textvariable=number2).grid(row=2, column=2)
entryNum3 = tk.Entry(root, textvariable=number3).grid(row=1, column=5)
entryNum4 = tk.Entry(root, textvariable=number4).grid(row=2, column=5)
call_result = partial(call_result, labelResult, labelTGsize, labelCGsize,labelpValue, number1, number2, number3, number4)
buttonCal = tk.Button(root, text="Calculate", command=call_result).grid(row=3, column=0)

new_mode = partial(new_mode, labelTitle, labelNum1, labelNum2, labelNum3, labelNum4, labelResult, labelTGsize, labelCGsize,labelpValue, number1, number2, number3, number4)
buttonMod = tk.Button(root, text="Change Mode", command=new_mode).grid(row=3, column=4)
root.mainloop()
