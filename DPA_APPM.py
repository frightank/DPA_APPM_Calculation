### Import Libraries ###
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

### Load Data from Text Files ###
# Define the path to the file
Pure_W_Rec_Ion_Filepath = r"C:\Users\j48032ja\OneDrive - The University of Manchester\Fusion CDT PhD\Project\Japan\Shimane University\SRIM Modelling\FC 5 keV 17 Degs\Pure W\E2RECOILIONIZ.csv"
Alloy_W_Rec_Ion_Filepath = r"C:\Users\j48032ja\OneDrive - The University of Manchester\Fusion CDT PhD\Project\Japan\Shimane University\SRIM Modelling\FC 5 keV 17 Degs\W-39Mo\E2RECOILIONIZ.csv"

# Load the data from the file
Pure_W_Rec_Ion_df = pd.read_csv(Pure_W_Rec_Ion_Filepath, encoding='unicode_escape', header=None, names=["DEPTH (Ang)", "Energy Absorbed by W (eV / Angstrom-Ion)", "Ioniastion by Recoils (eV / Angstrom-Ion)", "He Atoms ((Atoms/cm3) / (Atoms/cm2))"])
Alloy_W_Rec_Ion_df = pd.read_csv(Alloy_W_Rec_Ion_Filepath, encoding='unicode_escape', header=None, names=["DEPTH (Ang)", "Energy Absorbed by W (eV / Angstrom-Ion)", "Ioniastion by Recoils (eV / Angstrom-Ion)", "He Atoms ((Atoms/cm3) / (Atoms/cm2))"])

# Convert to lists
Depth_A_Pure_W = Pure_W_Rec_Ion_df["DEPTH (Ang)"].to_numpy()
Energy_Absorbed_by_Pure_W = Pure_W_Rec_Ion_df["Energy Absorbed by W (eV / Angstrom-Ion)"].to_numpy()
Ionisation_by_Recoils_Pure_W = Pure_W_Rec_Ion_df["Ioniastion by Recoils (eV / Angstrom-Ion)"].to_numpy()
He_Atoms_Pure_W = Pure_W_Rec_Ion_df["He Atoms ((Atoms/cm3) / (Atoms/cm2))"].to_numpy()

Depth_A_Alloy = Alloy_W_Rec_Ion_df["DEPTH (Ang)"].to_numpy()
Energy_Absorbed_by_Alloy = Alloy_W_Rec_Ion_df["Energy Absorbed by W (eV / Angstrom-Ion)"].to_numpy()
Ionisation_by_Recoils_Alloy = Alloy_W_Rec_Ion_df["Ioniastion by Recoils (eV / Angstrom-Ion)"].to_numpy()
He_Atoms_Alloy = Alloy_W_Rec_Ion_df["He Atoms ((Atoms/cm3) / (Atoms/cm2))"].to_numpy()

### Make adjustments to the data ###
# Convert the depth to nm
Depth_nm_Pure_W = Depth_A_Pure_W * 0.1
Depth_nm_Alloy = Depth_A_Alloy * 0.1

# Calculate Tdam
Tdam_Pure_W = (Energy_Absorbed_by_Pure_W - Ionisation_by_Recoils_Pure_W)
Tdam_Alloy = (Energy_Absorbed_by_Alloy - Ionisation_by_Recoils_Alloy)

# Material Constants
Avogadro = 6.022e23 # atoms / mol
Alloy_Conc_frac = 0.39 # at %
Atomic_Weight_W = 183.84 # g/mol
Atomic_Weight_Mo = 95.95 # g/mol
Molar_Volume_W = 9.55  # cm³/mol
Molar_Volume_Mo = 9.334  # cm³/mol
Energy_Disp_W = 90 # eV
Energy_Disp_Mo = 60 # eV
Atomic_Dens_W = 6.338e22 # atoms / cm^3
Atomic_Dens_Mo = 6.022e22 # atoms / cm^3

# Material Calculation
Atomic_Weight_Ave = (Atomic_Weight_W * (1 - Alloy_Conc_frac) + Atomic_Weight_Mo * Alloy_Conc_frac) # g/mol
Molar_Volume_Alloy = (Molar_Volume_W * (1 - Alloy_Conc_frac) + Molar_Volume_Mo * Alloy_Conc_frac) # cm³/mol
Energy_Disp_Alloy = (Energy_Disp_W * (1 - Alloy_Conc_frac) + Energy_Disp_Mo * Alloy_Conc_frac) # eV
Atomic_Dens_Alloy = (Avogadro / Molar_Volume_Alloy) # atoms / cm^3

# Beam Constants
joules_per_eV = 1.60218e-19 # J / eV
eV_per_J = 1 / joules_per_eV # eV / J
Ion_fluence_m2 = 1e20 # ions / m^2
ion_fluence_cm2 = Ion_fluence_m2 * 1e-4 # ions / cm^2

# Calculate vacancies
V_NRT_per_Ang_Ion_Pure_W = (0.8*Tdam_Pure_W) / (2*Energy_Disp_W)
V_NRT_per_cm_Ion_Pure_W = V_NRT_per_Ang_Ion_Pure_W * 1e8
V_NRT_per_Ang_Ion_Alloy = (0.8*Tdam_Alloy) / (2*Energy_Disp_Alloy)
V_NRT_per_cm_Ion_Alloy = V_NRT_per_Ang_Ion_Alloy * 1e8

# Calculate DPA
DPA_Pure_W = (V_NRT_per_cm_Ion_Pure_W * ion_fluence_cm2) / Atomic_Dens_W
DPA_Alloy = (V_NRT_per_cm_Ion_Alloy * ion_fluence_cm2) / Atomic_Dens_Alloy

# Convert He concentration to atoms/cm³
He_Atoms_per_cm3_Pure_W = He_Atoms_Pure_W * ion_fluence_cm2
He_Atoms_per_cm3_Alloy = He_Atoms_Alloy * ion_fluence_cm2

# Calculate appm of He
He_appm_Pure_W = (He_Atoms_per_cm3_Pure_W / Atomic_Dens_W) * 1e6
He_appm_Alloy = (He_Atoms_per_cm3_Alloy / Atomic_Dens_Alloy) * 1e6

# print(Energy_Disp_Alloy)
# print("Average Atomic Weight is", Atomic_Weight_Ave)
# print("Alloy Molar Volume is", Molar_Volume_Alloy)
# print("Alloy Atomic Density is", Atomic_Dens_Alloy)

# print(He_appm_Alloy)
# print(DPA_Alloy)



# Create a figure with two subplots, each of which has a left and right y-axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# First plot: Unalloyed W
color = 'tab:red'
ln1 = ax1.plot(Depth_nm_Pure_W, DPA_Pure_W, color=color, label="DPA")
ax1.set_xlabel("Specimen Depth from Irradiated surface (nm)", fontsize=12)
ax1.set_ylabel("Dose (DPA)", color=color, fontsize=12)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_title("Unalloyed W")
ax1.set_ylim(auto=True,top=0.13)

ax1_twin = ax1.twinx()
color = 'tab:blue'
ln2 = ax1_twin.plot(Depth_nm_Pure_W, He_appm_Pure_W, color=color, label="He appm", linestyle="--")
ax1_twin.set_ylabel("He Concentration (appm)", color=color, fontsize=12)
ax1_twin.tick_params(axis='y', labelcolor=color)
ax1_twin.set_ylim(auto=True, top=42500)

# Combine legends for both axes
lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='upper right')


# Second plot: W-39Mo Alloy
color = 'tab:red'
ax2.plot(Depth_nm_Alloy, DPA_Alloy, color=color, label="DPA")  # Replace Depth_nm and DPA with Depth_nm_Alloy and DPA_Alloy if different data
ax2.set_xlabel("Specimen Depth from Irradiated surface (nm)", fontsize=12)
ax2.set_ylabel("Dose (DPA)", color=color, fontsize=12)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_title("W-39at.%Mo Alloy")
ax2.set_ylim(auto=True,top=0.13)

ax2_twin = ax2.twinx()
color = 'tab:blue'
ax2_twin.plot(Depth_nm_Alloy, He_appm_Alloy, color=color, label="He appm", linestyle="--")  # Replace Depth_nm and He_appm with Depth_nm_Alloy and He_appm_Alloy if different data
ax2_twin.set_ylabel("He Concentration (appm)", color=color, fontsize=12)
ax2_twin.tick_params(axis='y', labelcolor=color)
ax2_twin.set_ylim(auto=True, top=42500)

# Adjust layout
fig.tight_layout()

# Show the plot
plt.show()

# # print max DPA and He appm
# print("Max DPA for Pure W is", max(DPA_Pure_W))
# print("Max DPA for W-39Mo Alloy is", max(DPA_Alloy))
# print("Max He appm for Pure W is", max(He_appm_Pure_W))
# print("Max He appm for W-39Mo Alloy is", max(He_appm_Alloy))

# # print depth of max DPA and He appm
# print("Depth of max DPA for Pure W is", Depth_nm_Pure_W[np.argmax(DPA_Pure_W)])
# print("Depth of max DPA for W-39Mo Alloy is", Depth_nm_Alloy[np.argmax(DPA_Alloy)])
# print("Depth of max He appm for Pure W is", Depth_nm_Pure_W[np.argmax(He_appm_Pure_W)])
# print("Depth of max He appm for W-39Mo Alloy is", Depth_nm_Alloy[np.argmax(He_appm_Alloy)])
