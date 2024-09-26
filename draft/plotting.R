library(rio)
library(dplyr)
library(ggplot2)
library(tidyr)
library(derivmkts)
options(scipen = 0)

iv_surface = rio::import("data/iv_surface_pred_price.csv") %>% mutate(erro_percentual = erro/price)
#iv_surface$r = 0
#iv_surface$d = 0


# GERANDO TABELAS ---------------------------------------------------------



iv_surface$type_option <- ifelse(iv_surface["moneyness"] < 0.99, "OTM",
                                 ifelse(iv_surface["moneyness"] > 1.01, "ITM", "ATM"))

iv_surface %>% group_by(type_option) %>% tally()

iv_surface <- iv_surface %>% mutate(media_erro = (price - pred_price)/pred_price)

tst <- iv_surface %>% group_by(moneyness, t_ano) %>% summarise(media_erro_money = mean(erro_percentual, na.rm = TRUE)) %>% round(digits = 2)
tst = tst %>% mutate(media_erro_money = paste0(media_erro_money*100, "%")) %>% spread(key = "moneyness", value = "media_erro_money") %>% mutate_all(~replace_na(., "Sem observações"))
write.csv(tst, file = "data/tabela_erro.csv", row.names = FALSE)

tst = iv_surface %>% group_by(type_option, t_ano) %>% tally()
tst = tst %>% spread(key = "type_option", value = "n")
tst = tst %>% round(digits = 2) %>% mutate(t_ano = t_ano)
  
###CALCULANDO VOLATILIDADE IMPLÍCITA
a <- lapply(1:nrow(iv_surface), function(n){
  
  iv_surface$r = 0
  iv_surface$d = 0
  
  
  iv_surface$iv_pred[n] <- derivmkts::bscallimpvol(price = iv_surface$price[n],
                          s = iv_surface$spot[n],
                          k = iv_surface$strike[n],
                          tt = iv_surface$t_ano[n],
                          r = iv_surface$r[n],
                          d = iv_surface$d[n])
  
  print(iv_surface$iv_pred[n])
  
})



# VISUALIZAÇÃO DAS CURVAS DE PREÇO ----------------------------------------



iv_surface = rio::import("data/iv_surface_pred_price.csv") %>%
  mutate(moneyness = spot / strike,
         t_ano = as.factor(round(t_ano, digits = 2))) #%>%
  #filter(t_ano > 0.5 & t_ano < 0.6)


ggplot(iv_surface %>%
         gather(key = "previsao_ou_mercado", value = "preco", -spot, -strike, -moneyness, -iv, -tipo, -t_ano, -erro) %>% 
         mutate(previsao_ou_mercado = ifelse(previsao_ou_mercado == "pred_price", "Preço do modelo", "Preço do mercado")),
       aes(x = moneyness, y = preco, shape = previsao_ou_mercado, color = previsao_ou_mercado)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('Preço do mercado' = "blue", 'Preço do modelo' = "red")) +
  scale_shape_manual(values = c('Preço do mercado' = 19, 'Preço do modelo' = 18)) +
  scale_alpha_manual(values = c("Preço do mercado" = 0.3, 'Preço do modelo' = 0.5)) +
  xlab("Moneyness") +
  ylab("Preço") +
  theme_bw() +
  theme(legend.title=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
        text = element_text(size=15),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom") +
  facet_grid(t_ano ~ .)


# VISUALIZAÇÃO DOS SMILES -------------------------------------------------



iv_surface = rio::import("data/iv_surface_pred_price.csv") %>%
  mutate(moneyness = spot / strike,
         t_ano = as.factor(round(t_ano, digits = 2))) #%>%
#filter(t_ano > 0.5 & t_ano < 0.6)


ggplot(iv_surface, aes(x = moneyness, y = iv)) +
  geom_point(size = 3) +
  xlab("Moneyness") +
  ylab("Volatilidade implícita") +
  theme_bw() +
  theme(legend.title=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
        text = element_text(size=15),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom") +
  facet_grid(t_ano ~ .)



# ggplot(iv_surface %>% gather(key = "previsao_ou_mercado", value = "preco", -spot, -strike, -moneyness, -iv, -tipo, -t_ano)) + +
#   geom_point(aes(x = moneyness, y = preco), size = 3) +
#   xlab("Moneyness") +
#   ylab("Volatilidade implícita") +
#   #ggtitle("Evolução da expectativa de inflação ao longo do ano de 2019") +
#   theme_bw() +
#   theme(legend.title=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
#         text = element_text(size=15),
#         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
