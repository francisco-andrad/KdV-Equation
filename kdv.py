import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Configurações ---
input_file = 'kdv_data.csv'
output_dir = 'images_kdv'

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
#         *** CALCULAR OS LIMITES GLOBAIS DOS EIXOS ***
# =============================================================================
# Encontra o valor mínimo e máximo de 'u' e 'x' em todo o arquivo
global_umin = df['u'].min()
global_umax = df['u'].max()
global_xmin = df['x'].min()
global_xmax = df['x'].max()

# Adiciona uma margem de 5% no eixo vertical para a onda não encostar no teto/chão do gráfico
margin_y = (global_umax - global_umin) * 0.05
if margin_y == 0: 
    margin_y = 0.1 # Segurança extra se a amplitude for exatamente 0

print(f"Escala vertical global definida: de {global_umin:.2f} a {global_umax:.2f}")
# =============================================================================

# 3. Obter a lista de todos os instantes de tempo únicos que foram salvos
times = df['time'].unique()
print(f"Encontrados {len(times)} instantes de tempo para plotar.")

# 4. Loop sobre cada instante de tempo para gerar um gráfico
for i, t in enumerate(times):
    print(f"Processando quadro {i+1}/{len(times)} (t = {t:.3f})...")
    
    # Isola apenas os dados do instante atual
    snapshot = df[df['time'] == t]
    
    # Ordena os valores por 'x' para garantir que a linha do plot não fique rabiscada
    # (útil caso o C++ tenha salvo fora de ordem por conta de concorrência)
    snapshot = snapshot.sort_values(by='x')
    
    plt.figure(figsize=(10, 6))
    
    # =========================================================================
    #   *** PLOT DA LINHA 1D ***
    # =========================================================================
    plt.plot(snapshot['x'], snapshot['u'], color='#1f77b4', linewidth=2)
    
    # Força os eixos a ficarem parados em todos os quadros
    plt.ylim(global_umin - margin_y, global_umax + margin_y)
    plt.xlim(-50, 50)
    # =========================================================================
    
    # Textos e estilo
    plt.title(f'Solução 1D no instante t = {t:.3f}')
    plt.xlabel('x')
    plt.ylabel('u(x, t)')
    
    # Um grid leve ajuda bastante a visualizar a velocidade da onda em vídeos
    plt.grid(True, linestyle='--', alpha=0.5) 
    
    output_path = os.path.join(output_dir, f'frame_{i:04d}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

print("\nProcesso concluído!")
print(f"Todos os {len(times)} quadros foram salvos no diretório '{output_dir}'.")
