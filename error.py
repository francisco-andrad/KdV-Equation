import pandas as pd
import matplotlib.pyplot as plt
import os

# --- Configurações ---
input_file = 'error_zk_data.txt'
output_dir = 'error_plots'

# --- Início do Script ---

# 1. Criar o diretório de saída se ele não existir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f"Diretório '{output_dir}' criado.")

# 2. Carregar os dados do arquivo de erro usando Pandas
try:
    print(f"Carregando dados do arquivo '{input_file}'...")
    df = pd.read_csv(input_file)
    print("Dados carregados com sucesso.")
except FileNotFoundError:
    print(f"ERRO: O arquivo '{input_file}' não foi encontrado. Verifique se a simulação C++ foi executada corretamente.")
    exit()

# --- 3. Gerar o Gráfico para o Erro da Massa ---
print("Plotando o erro da massa...")
# plt.figure(figsize=(10, 6))
plt.plot(df['time'], df['mass_error'])
# plt.title('Erro com relação à massa inicial')
plt.xlabel('t')
plt.ylabel('error')
plt.grid(True,alpha=0.6)
# plt.yscale('log') # Usar escala logarítmica é bom para ver pequenas variações
plt.legend()
plt.tight_layout()
# Salvar a imagem
mass_error_path = os.path.join(output_dir, 'mass_error.png')
plt.savefig(mass_error_path, dpi=150)
plt.close() # Fecha a figura para começar a próxima
print(f"Gráfico do erro da massa salvo em: '{mass_error_path}'")


# --- 4. Gerar o Gráfico para o Erro da Energia ---
print("Plotando o erro da energia...")
# plt.figure(figsize=(10, 6))
plt.plot(df['time'], df['energy_error'])
# plt.title('Erro com relação à energia inicial')
plt.xlabel('t')
plt.ylabel('error')
plt.grid(True,  alpha=0.6)
# plt.yscale('log') # Escala logarítmica também é útil aqui
plt.legend()
plt.tight_layout()
# Salvar a imagem
energy_error_path = os.path.join(output_dir, 'energy_error.png')
plt.savefig(energy_error_path, dpi=150)
plt.close()
print(f"Gráfico do erro da energia salvo em: '{energy_error_path}'")

print("\nProcesso concluído!")
