import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Configurações ---
input_file = 'zk_data.csv'
output_dir = 'images_zk'

# --- Início do Script ---

# 1. Criar o diretório de saída se ele não existir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Diretório '{output_dir}' criado.")

# 2. Carregar TODOS os dados da simulação de uma só vez usando Pandas
try:
    print(f"Carregando dados do arquivo '{input_file}'...")
    df = pd.read_csv(input_file)
    print("Dados carregados com sucesso.")
except FileNotFoundError:
    print(f"ERRO: O arquivo '{input_file}' não foi encontrado. Verifique se a simulação C++ foi executada corretamente.")
    exit()

# =============================================================================
#         *** MUDANÇA 1: CALCULAR OS LIMITES GLOBAIS AQUI ***
# =============================================================================
# Encontra o valor mínimo e máximo de 'u' em todo o arquivo de dados
global_vmin = df['u'].min()
global_vmax = df['u'].max()
print(f"Escala de cores global definida: de {global_vmin:.2f} a {global_vmax:.2f}")
# =============================================================================


# 3. Obter a lista de todos os instantes de tempo únicos que foram salvos
times = df['time'].unique()
print(f"Encontrados {len(times)} instantes de tempo para plotar.")


# 4. Loop sobre cada instante de tempo para gerar um gráfico
for i, t in enumerate(times):
    print(f"Processando quadro {i+1}/{len(times)} (t = {t:.3f})...")
    
    snapshot = df[df['time'] == t]
    grid = snapshot.pivot_table(index='y', columns='x', values='u')
    
    plt.figure(figsize=(10, 8))
    
    # =========================================================================
    #   *** MUDANÇA 2: USAR OS LIMITES GLOBAIS NO PLOT ***
    # =========================================================================
    # Adiciona os argumentos vmin e vmax à função imshow
    plt.imshow(grid, extent=[grid.columns.min(), grid.columns.max(), grid.index.min(), grid.index.max()],
               origin='lower', aspect='auto', cmap='viridis',
               vmin=global_vmin, vmax=global_vmax)
    # =========================================================================
    
    plt.colorbar(label='u(x,y;t)')
    # plt.title(f'Solução no instante t = {t:.3f}')
    plt.ylim(-20.0, 20.0)
    plt.xlim(-20.0,20.0)
    plt.xlabel('x')
    plt.ylabel('y')
    
    output_path = os.path.join(output_dir, f'frame_{i:04d}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

print("\nProcesso concluído!")
print(f"Todos os {len(times)} quadros foram salvos no diretório '{output_dir}'.")
