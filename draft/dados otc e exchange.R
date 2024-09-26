library(dplyr)
library(rio)
library(tidyverse)

otc <- rio::import("dados_tcc/WEBSTATS_OTC_DATAFLOW_csv_col.csv")
otc <- otc %>% filter(row_number() == 1) %>% select(30:70) %>% gather(key = "ano", value =  "valor")
otc <- otc %>% filter(grepl("-S2", ano))
otc$ano <- gsub("-S2", "", otc$ano)
otc <- otc %>% mutate(valor = as.numeric(valor)/1000000)
otc$tipo <- "Mercado de balcão"


exchange <- rio::import("dados_tcc/dados_bolsa_derivativos.csv")
exchange <- exchange %>%
gather(key = "ano", value = "valor") %>%
filter(grepl("-Q4", ano)) %>%
mutate(valor = as.numeric(valor)/1000000) %>% 
group_by(ano) %>%
summarise(valor = sum(valor))

exchange$ano <- gsub("-Q4", "", exchange$ano)
exchange <- exchange %>% filter(ano > "1997", ano < "2018")
exchange$tipo <- "Negociado em bolsa de valores"

total <- rbind(otc, exchange)

ggplot(total, aes(x = ano, y = valor, group = tipo)) +
  geom_line(aes(linetype = tipo), size = 1) +
  xlab("Ano") +
  ylab("Valor negociado (em trilhões US$)") +
  ggtitle("Volume negociado de derivativos no mercado de balcão") +
  theme_bw() +
  theme(legend.title=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

