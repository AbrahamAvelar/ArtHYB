library(dplyr)
library(tidyr)
library(ggplot2)

archivo="C:/Users/HP/Dropbox/Posdoc/Hybrids/Pct_of_hybrids/Sacc720.csv"
df=read.csv(archivo)

tousedf= df[!is.na(df$Pct_Hybrids) == 1,]
tousedf$GenomSpecies <- factor(tousedf$GenomSpecies, levels = c("SACE", "Hybrid", "SAPA"))

# FUNCTION
realizar_muestreo_y_graficar <- function(data, columna_separar, columna_agrupacion, columna_graficar, replicas) {
  # Crear todas las combinaciones posibles entre las categorías
  all_combinations <- expand.grid(
    separar = unique(data[[columna_separar]]),
    graficar = unique(data[[columna_graficar]])
  )
  names(all_combinations) <- c(columna_separar, columna_graficar)
  
  # Calcular el número de cuentas únicas de columna_agrupacion por cada columna_separar
  unique_counts <- data %>%
    group_by(.data[[columna_separar]]) %>%
    summarise(unique_count = n_distinct(.data[[columna_agrupacion]]), .groups = "drop")
  
  # Guardar resultados de las muestras aleatorias
  sampling_results <- list()
  
  for (i in 1:replicas) {
    # Muestreo aleatorio por la columna de agrupación
    random_sampled <- data %>%
      group_by(.data[[columna_agrupacion]]) %>%
      slice_sample(n = 1, replace = FALSE) %>%
      ungroup()
    
    # Calcular los conteos por las columnas especificadas
    sampled_counts <- random_sampled %>%
      group_by(.data[[columna_separar]], .data[[columna_graficar]]) %>%
      summarise(count = n(), .groups = "drop")
    
    # Unir con todas las combinaciones para agregar los ceros
    sampled_combined <- all_combinations %>%
      left_join(sampled_counts, by = c(columna_separar, columna_graficar)) %>%
      replace_na(list(count = 0))
    
    # Calcular los porcentajes
    sampled_combined <- sampled_combined %>%
      group_by(.data[[columna_separar]]) %>%
      mutate(percentage = count / sum(count) * 100) %>%
      ungroup()
    
    # Guardar resultados
    sampling_results[[i]] <- sampled_combined
  }
  
  # Calcular promedio y desviación estándar de los porcentajes
  sampling_avg_sd <- bind_rows(sampling_results) %>%
    group_by(.data[[columna_separar]], .data[[columna_graficar]]) %>%
    summarise(
      avg_percentage = mean(percentage, na.rm = TRUE),
      sd_percentage = sd(percentage, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Combinar los conteos únicos con los resultados finales
  sampling_avg_sd <- sampling_avg_sd %>%
    left_join(unique_counts, by = columna_separar)
  
  sampling_avg_sd <- sampling_avg_sd %>%
    mutate(
      x_label = paste0(.data[[columna_separar]], " (n=", unique_count, ")")
    )
  
  
    b<-  ggplot(sampling_avg_sd, aes(x = x_label, y = avg_percentage, fill = .data[[columna_graficar]])) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7) +
    geom_errorbar(aes(ymin = avg_percentage - sd_percentage, ymax = avg_percentage + sd_percentage),
                  width = 0.25, position = position_dodge(0.7)) +
    labs(
      x = columna_separar,
      y = "Average Percentage",
      title = paste("Average Percentage of", columna_graficar, "by", columna_separar,
                    "\n(Tomando 1 por", columna_agrupacion, ", réplicas =", replicas, ")")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    )
  plot(b)
  a<-ggplot(sampling_avg_sd, aes(x = "", y = avg_percentage, fill = .data[[columna_graficar]])) +
    geom_bar(stat = "identity", width = 1) +  # Crear las barras que serán convertidas en un pastel
    coord_polar(theta = "y") +  # Convertir a gráfico circular (pastel)
    labs(
      title = paste("Distribución de", columna_graficar, "por", columna_separar,
                    "\n(Tomando 1 por", columna_agrupacion, ", réplicas =", replicas, ")")
    ) +
    facet_wrap(~ .data[[columna_separar]], scales = "free") +  # Crear un gráfico de pastel por cada valor de columna_separar
    geom_text(aes(label = paste0(round(avg_percentage, 1), "%")),  # Añadir texto con el porcentaje
              position = position_stack(vjust = 0.5),  # Posicionar el texto en el centro de cada porción
              color = "white", size = 4) +  # Cambiar color del texto a blanco y ajustar tamaño
    theme_minimal() +
    theme(
      axis.text = element_blank(),  # Eliminar etiquetas del eje
      axis.title = element_blank(),  # Eliminar título del eje
      plot.title = element_text(hjust = 0.5),  # Centrar el título
      legend.title = element_blank()  # Eliminar título de la leyenda
    )
  plot(a)
  
}
#End FUNCTION


#Cuántas hay en total?
pdf("PCT_DistVSNat_Todas.pdf")
set.seed(2)
realizar_muestreo_y_graficar(  data = tousedf, columna_separar = "LocationType", 
                               columna_agrupacion = "Plate_Pos", columna_graficar = "GenomSpecies", replicas = 20) 
dev.off()

pdf("PCT_DistVSNat_Todas_IDFASTA.pdf")
set.seed(2)
realizar_muestreo_y_graficar(  data = tousedf, columna_separar = "LocationType", 
                               columna_agrupacion = "ID_Fasta", columna_graficar = "GenomSpecies", replicas = 20) 
dev.off()


#Si solo tomamos una por lccation
pdf("PCT_PorLocation_y_PorOBJID.pdf")
set.seed(2)
realizar_muestreo_y_graficar(  data = tousedf, columna_separar = "LocationType", 
  columna_agrupacion = "Location", columna_graficar = "GenomSpecies", replicas = 50) 
set.seed(2)
realizar_muestreo_y_graficar(  data = tousedf, columna_separar = "Region",columna_agrupacion = "Location",
                               columna_graficar = "GenomSpecies",replicas = 50) 
set.seed(2)
#Si solo tomamos una por OBJ_ID
realizar_muestreo_y_graficar(  data = tousedf,columna_separar = "LocationType",columna_agrupacion = "OBJ_ID",
                               columna_graficar = "GenomSpecies",replicas = 50)
set.seed(2)
realizar_muestreo_y_graficar(  data = tousedf,columna_separar = "Region",columna_agrupacion = "OBJ_ID",
                               columna_graficar = "GenomSpecies",replicas = 50)
dev.off()



subtouse=tousedf[ tousedf$LocationType=="Distillery",]
pdf("PCT_PorLocation_y_PorOBJID_OnlyDistillery.pdf")
set.seed(2)
realizar_muestreo_y_graficar( data = subtouse,columna_separar = "Region",columna_agrupacion = "Location",
                              columna_graficar = "GenomSpecies",replicas = 50)

realizar_muestreo_y_graficar( data = subtouse,columna_separar = "Region",columna_agrupacion = "OBJ_ID",
                              columna_graficar = "GenomSpecies",replicas = 50)
dev.off()
