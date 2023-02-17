# -*- coding: utf-8 -*-
import streamlit as st
import pandas as pd
import numpy as np
import math
import scipy.special as scsp
from functools import partial

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

def main():
	st.sidebar.title("What to do")
	app_mode = st.sidebar.selectbox("Choose the app mode", 
									["Significance calculator", "Minimum control group calculator"])
	if app_mode == "Significance calculator":
		significance_calculator()
	elif app_mode == "Minimum control group calculator":
		min_cg_calculator()


@st.cache
def z2p(z):
    """From z-score return p-value."""
    return 2*(1- (0.5 * (1 + scsp.erf(abs(z) / math.sqrt(2)))))

@st.cache
def zscore(p1, p2, n1, n2): # p1, p2 proportions
	p = (p1*float(n1) + p2*float(n2))/(float(n1) + float(n2))
	numerator = p1 - p2
	denominator = math.sqrt(p*(1-p)*((1/n1)+(1/n2)))
	return numerator/(denominator + 10**(-20))

def significance_calculator():
	st.title('Uplift & p-value calculator')
	tg_size = st.number_input('Target group counts', min_value = 1, max_value = 99999999, value = 200, key = '0101')
	tg_acceptors = st.number_input('Target group acceptors', min_value = 0, max_value = 99999999, value = 0, key ='0201')
	cg_size = st.number_input('Control group counts', min_value = 1, max_value = 99999999, value = 100, key ='0301')
	cg_acceptors = st.number_input('Control group acceptors', min_value = 0, max_value = 99999999, value = 0, key ='0401')

	z_score = zscore(tg_acceptors/float(tg_size),cg_acceptors/float(cg_size),float(tg_size),float(cg_size))
	st.write('Uplift is ', round(float(tg_acceptors)/tg_size - float(cg_acceptors)/cg_size,2)*100)
	st.write('Z-score is ', round(z_score, 2))
	st.write('p-value is ', round(z2p(z_score),2))

def min_cg_calculator():
	st.title('Minimum control group calculator')

	epsilon = 0.001
	step = 0.005
	initial_solution = 0.1 # Initial split (for CG from total)
	solution = []

	total_counts = st.number_input('Total counts', min_value = 1, max_value = 99999999, value = 200, key = '0501')
	expected_uplift = st.number_input('Expected uplift (%)', min_value = 0.00, max_value = 100.0, value = 1.0, step = 0.05, key ='0601')/100
	significance = st.number_input('Significance (%)', min_value = 80.0, max_value = 100.0, value = 90.0, step = 0.1, key ='0701')/100
	max_acceptance_cg = st.number_input('Maximum acceptance control group (%)', min_value = 0.0, max_value = 100.0, value = 1.0, step = 0.05, key ='0801')/100

	def evaluate(value):
		cg_acceptance = max_acceptance_cg
		tg_acceptance = expected_uplift + cg_acceptance
		tgsize = (1-value)*total_counts
		cgsize = value*total_counts

		if (z2p(zscore(tg_acceptance, cg_acceptance, tgsize, cgsize )) >= 1 - significance + epsilon): #and cg_acceptors < 1
			return False
		if max_acceptance_cg*value*total_counts < 2: #minimum number of records
			return False
		else:
			return True

	def descend(init_value):
		value = init_value - step
		solution.append(value)
		if evaluate(value) and value > 0.02:
			print("Descending to find a better solution")
			descend(value)
		elif (value <= 0.02) or (not evaluate(value)):
			solution.append(value)

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
			cg_size = int(solution[len(solution)-1]*total_counts)
			tg_size = total_counts-cg_size
			print((solution[0],tg_size,cg_size))
			return (1,solution[0],tg_size,cg_size)
		else:
			ascend(initial_solution)
			if len(solution) > 0:
				#Show TG/CG size, TG/CG acceptors and p-value
				cg_size = round(solution[len(solution)-1]*total_counts)
				tg_size = total_counts-cg_size
				print((solution[0],tg_size,cg_size))
				return (1,solution[0],tg_size,cg_size)

			else:
				print("No solution")
				return (0,0, total_counts, 0)

	if calculate()[0] == 1:
		st.write('Final solution: ', calculate()[1])
		st.write('Maxium Target group size: ', calculate()[2])
		st.write('Minimum Control group size: ', calculate()[3])
	else:
		st.write('No solution')


if __name__ == "__main__":
    main()