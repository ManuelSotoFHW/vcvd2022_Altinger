import cmath
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

# Constants
g = 9.81
mu = 0.5

# Coefficients of Friction According to Road Types and Road Conditions
Coefficients_Of_Friction = pd.DataFrame([[0.5, 0.35, 0], [0.15, 0.08, 0], [0, 0, 0.05], [0.35, 0, 0], [0.3, 0, 0]],
                                        index=('concrete', 'ice', 'water', 'gravel', 'sand'),
                                        columns=('dry', 'wet', 'aquaplaning'))

print(Coefficients_Of_Friction)


def Break_Distance(V, ag, m, angle):
    return pow(V, 2) / (2 * ag * (m * math.cos(angle) + math.sin(angle)))


def Accelaration(V, s):
    return -pow(V, 2) / (2 * s)


def Break_time(V, aa, s):
    a = 0.5 * aa
    b = V
    c = -s
    dis = (b ** 2) - (4 * a * c)
    t = (-b - cmath.sqrt(dis)) / (2 * a)
    # t2 = (-b + cmath.sqrt(dis)) / (2 * a)
    return abs(t)


def Rule_Of_Thumb_Normal(V):
    return pow((V / 10), 2)


def Rule_Of_Thumb_Danger(V):
    return Rule_Of_Thumb_Normal(V) * (1 / 3)


def Rule_Of_Thumb_Reaction(V):
    return (V / 10) * 3


V0 = float(input("Please enter the initial velocity in [m/s]: "))  # 14
alpha = float(input("Please enter the road inclination in degrees: "))
print(" ")
print("Please select one of Road Type options: concrete, ice, water, gravel, sand ")
Road_Type = str.casefold(input("Please type your Road Type selection: "))
print(" ")
print("Please select one of Road Conditions options: dry, wet, aquaplaning ")
Road_Condition = str.casefold(input("Please type your Road Condition selection: "))

mu = float(Coefficients_Of_Friction.loc[Road_Type][Road_Condition])
print(" ")
print("The dynamic coefficient friction for the selected road type and condition is mu=", mu)

# Convert road angle in degrees to radians
alpha = alpha * (math.pi / 180)

# Calculation od the acceleration
acc = -g * (mu * math.cos(alpha) + math.sin(alpha))

# Calculation of Breaking distance for simply math model
s_simply_math_model = Break_Distance(V0, g, mu, alpha)
print("The breaking distance is s=", s_simply_math_model, " [m]")

# Calculation of Breaking time for simply math model
t_simply_math_model = Break_time(V0, acc, s_simply_math_model)
print("The breaking time is t= ", t_simply_math_model, " [s]")

# "Sampling" of the breaking time for simply math model
Time_simply_math_model = np.arange(0, t_simply_math_model, t_simply_math_model / 2000)

# Calculation of the velocity in each time step for simply math model
Velocity_simply_math_model = V0 + acc * Time_simply_math_model

# Calculation of the distance in each time step for simply math model
Distance_simply_math_model = 0 + (V0 * Time_simply_math_model) + (0.5 * acc * np.square(Time_simply_math_model))

# Plots
fig, axs = plt.subplots(2)
fig.set_size_inches(8.3, 11.7)
fig.suptitle("Breaking Simulation")

axs[0].plot(Time_simply_math_model, Velocity_simply_math_model, 'tab:red')
axs[0].set_xlabel("Time [s]")
axs[0].set_ylabel("Velocity [m/s]")
axs[0].grid()
# axs[0].show()

axs[1].plot(Time_simply_math_model, Distance_simply_math_model, 'tab:green')
axs[1].set_xlabel("Time [s]")
axs[1].set_ylabel("Distance [m]")
axs[1].grid()

plt.savefig("Breaking_Simulation.pdf", format="pdf", bbox_inches="tight")
plt.show()
