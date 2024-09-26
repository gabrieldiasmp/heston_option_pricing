import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D  
from matplotlib import cm

# iv_surface = pd.read_csv("data/spy_options.csv")
# iv_surface.rename(columns={'Strike': 'strike',
#                            'Implied Volatility': 'iv',
#                            'Last Price': "price"}, inplace=True)

# iv_surface["moneyness"] = iv_surface["spot"] / iv_surface["strike"]

# iv_surface = iv_surface.loc[(iv_surface["tipo"] == "calls") & 
# (iv_surface["t_ano"] > 0.1) & 
# (iv_surface["t_ano"] < 5) & 
# (iv_surface["moneyness"] < 1.66) & 
# (iv_surface["moneyness"] > 0.7) & 
# (iv_surface["iv"] > 0) &
# (iv_surface["strike"].isin([180, 200, 220, 240, 270, 299, 315, 330, 350, 370, 400, 420])), 
# ["spot", "strike", 'moneyness', "price", "iv", "tipo", "t_ano"]]


iv_surface = pd.read_csv("data/iv_surface_pred_price.csv")

###PLOT DOS PONTOS
fig = plt.figure()
plt.style.use('default')
fig.set_facecolor("w")
ax = plt.axes(projection='3d')
ax.scatter3D(iv_surface["moneyness"], iv_surface["t_ano"], iv_surface["price"], label = "Preço de mercado", alpha = 0.8)
ax.scatter3D(iv_surface["moneyness"], iv_surface["t_ano"], iv_surface["pred_price"], c = "red", label = 'Preço do modelo', alpha = 0.3, marker = "s")
ax.view_init(elev=30, azim = 115) #azim=45
ax.set_xlabel('Moneyness')
ax.set_ylabel('Prazo de expiração')
ax.set_zlabel('Preço da opção')
ax.legend()

###PLOT DA SUPERFÍCIE DOS PREÇOS PREDITOS
fig = plt.figure()
plt.style.use('default')
fig.set_facecolor("w")
ax = Axes3D(fig)
surf = ax.plot_trisurf(iv_surface["moneyness"], iv_surface["t_ano"], iv_surface["pred_price"], cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('Moneyness')
ax.set_ylabel('Prazo de expiração')
ax.set_zlabel('Preço')
ax.view_init(elev=30., azim=115)
plt.show()

###PLOT DA SUPERFÍCIE DA VOLATILIDADE
fig = plt.figure()
plt.style.use('default')
fig.set_facecolor("w")
ax = Axes3D(fig)
surf = ax.plot_trisurf(iv_surface["moneyness"], iv_surface["t_ano"], iv_surface["iv"], cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('Moneyness')
ax.set_ylabel('Prazo de expiração')
ax.set_zlabel('Volatilidade implícita')
ax.view_init(elev=30., azim=115)
plt.show()
