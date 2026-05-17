import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Configurações ---
input_file = 'centers_data.csv'
output_dir = 'monotonicity_plots'

# --- Início do Script ---

# 1. Criar o diretório de saída se ele não existir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Diretório '{output_dir}' criado.")

# 2. Carregar os dados do arquivo de centros
try:
    print(f"Carregando dados do arquivo '{input_file}'...")
    df = pd.read_csv(input_file)
    # Renomeia as colunas para facilitar o acesso
    df.columns = ['time', 'sx', 'sy']
    print("Dados carregados com sucesso.")
except FileNotFoundError:
    print(f"ERRO: O arquivo '{input_file}' não foi encontrado. Verifique se a simulação C++ foi executada corretamente.")
    exit()
except ValueError:
     print(f"ERRO: O arquivo '{input_file}' parece ter um número incorreto de colunas. Verifique o cabeçalho no código C++.")
     exit()

# --- 3. Gráfico 1: Plotar as diferenças dos centros ---
print("Gerando gráfico das componentes de separação...")
# plt.figure(figsize=(10, 6))
plt.plot(df['time'], df['sx'], label=' $\\langle x \\rangle_E - \\langle x \\rangle_M$')
plt.plot(df['time'], df['sy'], label='$ \\langle y \\rangle_E - \\langle y \\rangle_M$')
# plt.title('Diferença entre a posição dos centros de Energia e Massa')
plt.xlabel('t')
plt.ylabel('difference')
plt.legend()
plt.grid(True, alpha=0.6)
plt.tight_layout()

# Salvar a primeira imagem
separation_path = os.path.join(output_dir, 'zk_mon_form.png')
plt.savefig(separation_path, dpi=150)
plt.close()
print(f"Gráfico salvo em: '{separation_path}'")


# --- 4. Gráfico 2: Calcular e plotar as derivadas ---
print("Calculando e plotando as derivadas...")

# Converter as colunas para arrays NumPy para facilitar o cálculo
t = df['time'].to_numpy()
sx = df['sx'].to_numpy()
sy = df['sy'].to_numpy()

# Calcular a derivada por diferenças finitas centradas
# f'(i) = (f(i+1) - f(i-1)) / (t(i+1) - t(i-1))
# Descartamos o primeiro e o último ponto, onde a fórmula não se aplica.
t_deriv = t[1:-1] # Novo eixo de tempo para as derivadas
dsx_dt = (sx[2:] - sx[:-2]) / (t[2:] - t[:-2])
dsy_dt = (sy[2:] - sy[:-2]) / (t[2:] - t[:-2])

# Criar o segundo gráfico
# plt.figure(figsize=(10, 6))
plt.plot(t_deriv, dsx_dt, label='x direction')
plt.plot(t_deriv, dsy_dt, label='y direction')
# plt.title('Derivada da fórmula de monotonicidade:')
plt.xlabel('t')
plt.ylabel('Derivative')
plt.legend()
plt.grid(True,alpha=0.6)
plt.tight_layout()

# Salvar a segunda imagem
derivative_path = os.path.join(output_dir, 'zk_mon_form_derv.png')
plt.savefig(derivative_path, dpi=150)
plt.close()
print(f"Gráfico salvo em: '{derivative_path}'")

print("\nProcesso concluído!")
