import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

EXPERIMENTS = 10

def create_heatmap():
    for i in range(EXPERIMENTS):
        file = np.loadtxt(f'./resultados/hipertermia_heterogeneo_np_{i+1}.txt')
        # Define the plot
        fig, ax = plt.subplots(figsize=(13,7))
        
        # Add title to the Heat map
        title = "Heat Map"
        
        # Set the font size and the distance of the title from the plot
        plt.title(title,fontsize=18)
        ttl = ax.title
        ttl.set_position([0.5,1.05])
        
        # Hide ticks for X & Y axis
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Remove the axes
        ax.axis('off')
        
        sns.heatmap(file, cmap='coolwarm', ax=ax, vmin=37, vmax=42)
        
        # Display the Pharma Sector Heatmap
        plt.savefig(f'./resultados/heatmap_{i+1}.png')
        
def plot_temp():
    muscle = np.zeros(EXPERIMENTS)
    tumor = np.zeros(EXPERIMENTS)
    sdv_muscle = np.zeros(EXPERIMENTS)
    sdv_tumor = np.zeros(EXPERIMENTS)
    skin = np.linspace(0, EXPERIMENTS, EXPERIMENTS)*10
    for i in range(EXPERIMENTS):
        file = np.loadtxt(f'./resultados/hipertermia_heterogeneo_np_{i+1}.txt')
        sdv = np.loadtxt(f'./resultados/hipertermia_sdv_{i+1}.txt')
        sdv_muscle[i] = sdv[25][25]
        sdv_tumor[i] = sdv[50][50]
        muscle[i] = file[25][25]
        tumor[i] = file[50][50]
    
    plt.figure()
    plt.plot(skin, muscle, '-k')
    plt.plot(skin, tumor, '-c')
    plt.fill_between(skin, (tumor + sdv_tumor), (tumor - sdv_tumor), alpha=0.5)
    plt.fill_between(skin, (muscle + sdv_muscle), (muscle - sdv_muscle), alpha=0.4)
    plt.xlabel('Tempo (min)')
    plt.ylabel('Temperatura (°C)')
    plt.legend(['Média da temperatura no tecido saudável', 'Média da temperatura no tumor'], loc=2, prop={'size': 8})
    plt.savefig('./resultados/temperaturas_medias.png')

create_heatmap()