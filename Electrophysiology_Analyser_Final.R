Hello

list.of.packages <- c("shiny","ggplot2","dplyr","ggpubr","shinythemes","shinyFiles",
                      "data.table","reticulate","tidyr","rstatix","factoextra","stats","cluster","shinycssloaders")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)
lapply(list.of.packages, require, character.only = TRUE)


# Define UI
ui <- fluidPage(
  theme = shinytheme("sandstone"),
  sidebarLayout(
    sidebarPanel(
      textInput("folder", "Enter Folder Path:", placeholder = "/path/to/folder"),
      fileInput("trigger_file", "Upload Trigger Times:", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      radioButtons("TriggerUnits", "Trigger Time Units:", 
                   choices = c("Seconds", "Samples"),
                   selected = "Seconds"),  # Default selection
      textInput("viewing_window", "Viewing Window:", value = "-0.05,0.05"),
      textInput("response_window_before_trigger", "Spontaneous Window:", value = "1,2"),
      textInput("response_window_after_trigger", "Response Window:", value = "0.003,0.02"),
      textInput("bin_width", "Bin Size:", value = "0.001"),
      selectInput("state", "Filtering:",
                  choices = c("No Filtering","Two Sided","Increased Activity","Decreased Activity")
      ),
      
      actionButton("generate_plot", "Generate plots"),
      actionButton("save_csv", "Save data (csv)"),
      actionButton("save_plot", "Save plot (pdf)"),
      width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("PSTH",
                 shinycssloaders::withSpinner(
                   plotOutput("PSTH"),type=1,color="blue")
        ),
        tabPanel("Raster Plots",
                 plotOutput("Raster_Plots")
        )
      ),
      verbatimTextOutput("dir")
    )
  )
)
# Define server logic
server <- function(input, output, session) {
  
  observeEvent(input$generate_plot, {


      # Server logic starts here after setting the working directory
      
      # Reactive expression to store the plot and data frames
      plot_data <- reactive({
        # Get input values
        folder_path <- input$folder
        setwd(folder_path)

        viewing_window <- as.numeric(unlist(strsplit(input$viewing_window, ",")))
        response_window_before_trigger <- as.numeric(unlist(strsplit(input$response_window_before_trigger, ",")))
        response_window_after_trigger <- as.numeric(unlist(strsplit(input$response_window_after_trigger, ",")))
        bin_width <- as.numeric(input$bin_width)
        filtering <- ifelse(input$state=="No Filtering", "no", "yes")
        
        spont_dur<-response_window_before_trigger[2]-response_window_before_trigger[1]
        evoked_dur<-response_window_after_trigger[2]-response_window_after_trigger[1]
        
        # Check if trigger file is uploaded
        
        # Read trigger file

        
        if (input$TriggerUnits == "Seconds") {
          trigger_times <- as.numeric(read.csv(input$trigger_file$datapath, sep = "", header = FALSE)[,1])
        } else {
          trigger_times <- as.numeric(read.csv(input$trigger_file$datapath, sep = "", header = FALSE)[,1]) / 3.003003003003003e+04
        }
        print(trigger_times)
        np <- import("numpy")
        depth<-as.data.frame(read.table("cluster_info.tsv",sep="\t",header=TRUE))
        depth$id<-depth$id+1
        depth$depth<-800-depth$depth
        
        spike_clusters <- np$load(paste0(folder_path,"/spike_clusters.npy"))+1
        spike_times <- (np$load(paste0(folder_path,"/spike_times.npy"))/3.003003003003003e+04)

        spike_times_and_clusters<-data.frame(spike_times,spike_clusters)
        cluster_id<-read.table(file = paste0(folder_path,'/cluster_group.tsv'), sep = '\t', header = TRUE)
        spike_times_and_clusters$spike_clusters<-spike_times_and_clusters$spike_clusters
        cluster_id$cluster_id<-cluster_id$cluster_id+1
        spike_times_and_clusters$group <- cluster_id$group[match(spike_times_and_clusters$spike_clusters, 
                                                                 cluster_id$cluster_id)]
        spike_times_and_clusters<-spike_times_and_clusters[spike_times_and_clusters$group!="noise",]
        
        
        clusters_all_spike_times <- list()
        for (neuron_id in unique(spike_times_and_clusters$spike_clusters)) {
          clusters_all_spike_times[[as.character(neuron_id)]] <- spike_times_and_clusters$spike_times[spike_times_and_clusters$spike_clusters == neuron_id]
        }
        
        clusters_within_window_spike_times <- list()
        for (j in 1:length(clusters_all_spike_times)) {
          spikes_within_window <- list()
          stim_numbers <- list()
          neuron_spikes <- clusters_all_spike_times[[j]]
          for (i in 1:length(trigger_times)) {
            within_window <- neuron_spikes[neuron_spikes > (trigger_times[i] - 11) &
                                             neuron_spikes < (trigger_times[i] + 11)] - trigger_times[i]
            if (length(within_window) == 0) {
              within_window <- NA  # Placeholder for no spikes
            }
            spikes_within_window[[i]] <- within_window
            stim_numbers[[i]] <- rep(i, length(within_window))
          }
          neuron_df <- data.frame(
            Cluster_ID = rep(names(clusters_all_spike_times)[j], length(unlist(spikes_within_window))), # Add ID column with neuron ID
            Spike_Time = unlist(spikes_within_window),
            Trial_Number = unlist(stim_numbers)
          )
          clusters_within_window_spike_times[[j]] <- neuron_df
        }
        
        clusters_within_window_spike_times <- clusters_within_window_spike_times[sapply(clusters_within_window_spike_times, nrow) > 0]
        
        cluster_info_per_trial_within_viewing_window <- list()
        for (j in 1:length(clusters_within_window_spike_times)){
          z <- data.frame(matrix(0, nrow = length(trigger_times), ncol = 3)) # Initialize dataframe z with an additional column
          
          for (i in 1:length(trigger_times)){
            z[i,1] <- sum(clusters_within_window_spike_times[[j]]$Spike_Time <= (response_window_before_trigger[1]*-1) & 
                            clusters_within_window_spike_times[[j]]$Spike_Time >= (response_window_before_trigger[2]*-1)&
                            clusters_within_window_spike_times[[j]]$Trial_Number == i,na.rm = TRUE)
            
            z[i,2] <- sum(clusters_within_window_spike_times[[j]]$Spike_Time <= (response_window_after_trigger[2]) & 
                            clusters_within_window_spike_times[[j]]$Spike_Time >= (response_window_after_trigger[1])&
                            clusters_within_window_spike_times[[j]]$Trial_Number == i,na.rm = TRUE)
          }
          
          colnames(z) <- c("Spikes_Before", "Spikes_After", "Cluster_ID")
          
          z$Cluster_ID <- unique(clusters_within_window_spike_times[[j]]$Cluster_ID)
          z$Cluster_ID<-as.numeric(z$Cluster_ID)
          z$Trial_Number<- 1:length(trigger_times)
          z$Rate_Before<-z$Spikes_Before/(response_window_before_trigger[2]-response_window_before_trigger[1])
          z$Rate_After<-z$Spikes_After/(response_window_after_trigger[2]-response_window_after_trigger[1])
          
          cluster_info_per_trial_within_viewing_window[[j]] <- z
        }
        
        column_means_list <- lapply(cluster_info_per_trial_within_viewing_window, function(df) {
          numeric_df <- df %>% select(-4)  # Exclude column 4
          colMeans(numeric_df, na.rm = TRUE)
        })
        cluster_info_averaged_trials <- do.call(rbind, column_means_list)
        cluster_info_averaged_trials<-as.data.frame(cluster_info_averaged_trials)
        
        
        if(input$state=="No Filtering"){
          alternative<-NULL
        }else if(input$state=="Two Sided"){
          alternative<-"two.sided"
        }else if(input$state=="Increased Activity"){
          alternative<-"greater"
        }else if(input$state=="Decreased Activity"){
          alternative<-"less"}
        
        statistic<-data.frame("Cluster_ID"=NA,"p_value"=NA)
        for (i in 1:nrow(cluster_info_averaged_trials)){
          compiled_long <- as.data.frame(cluster_info_per_trial_within_viewing_window[i]) %>%
            gather(key = "trigger", value = "Rate", Rate_Before, Rate_After )
          statistic[i,1]<-compiled_long$Cluster_ID[1]
          statistic[i,2]<-wilcox_test(data=compiled_long,Rate ~ trigger, paired = TRUE, alternative = alternative)$p
        }
        
        significant_clusters<-statistic[statistic$p_value<0.05,1]
        significant_clusters <- as.numeric(significant_clusters[!is.na(significant_clusters)])
        
        cluster_info_averaged_trials <- cluster_info_averaged_trials %>%
          mutate(Significant = if_else(Cluster_ID %in% significant_clusters, TRUE, FALSE))
        cluster_info_averaged_trials$Increased_Activity<-cluster_info_averaged_trials$Rate_After>cluster_info_averaged_trials$Rate_Before
        
        cluster_info_averaged_trials_significant<-cluster_info_averaged_trials[cluster_info_averaged_trials$Significant==TRUE,]
        
        significant_clusters_within_window_spike_times <- lapply(clusters_within_window_spike_times, function(df) {
          if (df$Cluster_ID[1] %in% significant_clusters) {
            return(df)
          } else {
            return(NULL)
          }
        })
        significant_clusters_within_window_spike_times <- Filter(Negate(is.null), significant_clusters_within_window_spike_times)
        
        
        filtered_cluster_info_averaged_trials_per_trial_within_viewing_window <- lapply(cluster_info_per_trial_within_viewing_window, function(df) {
          if (df$Cluster_ID[1] %in% significant_clusters) {
            return(df)
          } else {
            return(NULL)
          }
        })
        filtered_cluster_info_averaged_trials_per_trial_within_viewing_window <- Filter(Negate(is.null), filtered_cluster_info_averaged_trials_per_trial_within_viewing_window)
        
        if (filtering == "yes") {
          Data <- significant_clusters_within_window_spike_times
        } else {
          Data <- clusters_within_window_spike_times
        }
        
        firing_rates <- list()
        for (j in 1:length(Data)) {
          neuron_df <- Data[[j]]
          
          bin_start_first <- viewing_window[1]
          bin_start_last <- viewing_window[2] - bin_width
          
          num_bins <- ceiling((bin_start_last - bin_start_first) / bin_width)
          
          spike_counts <- numeric(num_bins)  # Initialize with zeros, not NAs
          
          for (i in 1:num_bins) {
            bin_start <- bin_start_first + (i - 1) * bin_width
            bin_end <- bin_start + bin_width
            spike_counts[i] <- sum(neuron_df$Spike_Time > bin_start & neuron_df$Spike_Time <= bin_end, na.rm = TRUE)
          }
          
          # Prevent division by zero
          max_spikes <- max(spike_counts, na.rm = TRUE)
          if (max_spikes == 0) {
            firing_rate <- rep(0, num_bins)  # Avoid division by zero
          } else {
            firing_rate <- spike_counts / max_spikes
          }
          
          bin_starts <- seq(bin_start_first, by = bin_width, length.out = num_bins)
          
          neuron_id <- rep(neuron_df[1,1], num_bins)
          bin_number <- 1:num_bins
          
          neuron_firing_df <- data.frame(
            Neuron_ID = neuron_id,
            Bin_Number = bin_number,
            Spike_Count = spike_counts,
            Firing_Rate = firing_rate,
            Bin_Start_Time = bin_starts
          )
          
          firing_rates[[j]] <- neuron_firing_df
        }
        
        all_firing_rates_df <- do.call(rbind, firing_rates)
        all_firing_rates_df$depth<-depth$depth[match(all_firing_rates_df$Neuron_ID,
                                                     depth$id)]
        all_firing_rates_df <- all_firing_rates_df[order(all_firing_rates_df$depth,decreasing = TRUE), ]
        all_firing_rates_df$DepthOrder<-rep(1:length(unique(all_firing_rates_df$Neuron_ID)),
                                            each=length(unique(all_firing_rates_df$Bin_Number)))
        
        
        viewing_window<-as.numeric(viewing_window)
        
        total_spike_counts <- rep(0, ceiling((viewing_window[2] - viewing_window[1]) / bin_width))
        bin_starts <- seq(viewing_window[1], by = bin_width, length.out = length(total_spike_counts))
        bin_numbers <- seq_along(total_spike_counts)
        
        for (neuron_df in Data) {
          spike_counts <- table(cut(as.numeric(neuron_df$Spike_Time), breaks = seq(viewing_window[1], viewing_window[2], bin_width)))
          total_spike_counts <- total_spike_counts + as.vector(spike_counts)
        }
        
        total_spike_counts_df <- data.frame(
          Bin_Number = bin_numbers,
          Bin_Start_Time = bin_starts,
          Total_Spike_Count = total_spike_counts
        )

        r<-data.frame()
        h<-data.frame()
        
        for (i in 1:length(Data)){
          r<-data.frame(Data[i])[,2:3]
          colnames(r)<-c("Spike_Time","Trial_Number")
          for (j in 1:length(trigger_times)){
            Data1<-r[r$Trial_Number==j,1]
            h[j,i]<-((sum(Data1>0)/evoked_dur)-(sum(Data1<0)/spont_dur))/((sum(Data1>0)/evoked_dur)+(sum(Data1<0)/spont_dur))
          }}
        
       
        SpikeRatio<-data.frame("SpikeRatio"=rowMeans(h,na.rm = TRUE),"Trial"=1:length(trigger_times))
        
        max_spike_count <- max(total_spike_counts_df$Total_Spike_Count)
        y_limit <- (max_spike_count * 1.2)
        
        k <- data.frame(
          group = c("Not modulated", "Modulated"),
          value = c(length(clusters_within_window_spike_times) - length(significant_clusters_within_window_spike_times),
                    length(significant_clusters_within_window_spike_times))
        )
        k <- k %>%
          arrange(desc(group)) %>%
          mutate(prop = value / sum(value) * 100) %>%
          mutate(ypos = cumsum(prop) - 0.5 * prop)
        
        
        Clustering_ZScore<-c(as.numeric((cluster_info_averaged_trials$Spikes_After-cluster_info_averaged_trials$Spikes_Before)/sd(cluster_info_averaged_trials$Spikes_Before)))
        silhouette_score<-function(l){
          km<-kmeans(Clustering_ZScore, centers=l, nstart = 25)
          ss<-silhouette(km$cluster,dist(Clustering_ZScore))
          mean(ss[,3])
        }
        l<-2:15
        avg_sil<-sapply(l,silhouette_score)
        silhouette_plot<-data.frame("Number"=1:length(avg_sil),
                                    "SilhouetteScore"=avg_sil)
        
        

        #PSTH, population PSTH, raster plots, and modulation index per trial

        col_names <- c("Cluster","SpikeTime", "Trial")
        window <- c(viewing_window[1], viewing_window[2])
          create_plot <- function(data) {
          Data <- data
          colnames(Data) <- col_names
          ggplot(Data, aes(x = SpikeTime, y = Trial)) +
            scale_shape_identity() +
            geom_point(shape = 20, size = 2, color = "red") +
            theme_classic() +
            ylim(c(-3, length(trigger_times)+3)) +
            xlim(window)+
            ggtitle(paste0("Unit ",data[1,1]))+
            theme(
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(color="darkgreen", size=9, face="bold.italic")) +
            ylab("Trial") +
            xlab("Time (s)") +
            geom_vline(xintercept = 0, color = "blue")
        }
        
        # Apply the function to each element of the list
          rasters <- lapply(significant_clusters_within_window_spike_times, create_plot)
          rasterplots <- ggarrange(plotlist = rasters)

          
          p1<-ggplot(data=silhouette_plot,aes(x=Number, y=SilhouetteScore))+
            geom_line()+
            geom_point()+
            theme_bw()+
            theme(panel.background = element_blank(),
                  plot.background = element_blank())+
            theme(text=element_text(size=21))
          
          
        p2<-ggplot(all_firing_rates_df, aes(x = Bin_Start_Time, y = DepthOrder, fill = Firing_Rate)) +
          geom_raster(interpolate = FALSE) +
          scale_fill_viridis_c(option = "B", na.value = "black") +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),  # Remove major gridlines
                panel.grid.minor = element_blank(),  # Remove minor gridlines
                panel.background = element_blank(),  # Remove panel background
                plot.background = element_blank(),   # Remove plot background
                axis.line = element_line(color = "black")) +  # Remove axis lines
          ylab("Neuron") +  # Corrected y-axis label syntax
          xlab("Time since stimulus (s)") +
          scale_x_continuous(breaks=seq(viewing_window[1],
                                        viewing_window[2],by=0.1))+
          geom_vline(xintercept = 0, color = "red")+
          theme(text=element_text(size=21))
        
        
        p3<-ggplot(total_spike_counts_df, aes(x = Bin_Start_Time, y = (Total_Spike_Count))) +
          geom_line() +
          geom_area(fill="blue") +  
          ylab("Total Spike Count") +
          ylim(c(0,300)) +
          xlab("") +
          theme_bw()+
          # geom_vline(xintercept = 5, color = "red", size = 1)+
          xlab("Time since stimulus (s)")+
          scale_x_continuous(breaks=seq(viewing_window[1],
                                        viewing_window[2],by=0.1))+
          geom_vline(xintercept = 0, color = "red")+
          theme(                panel.background = element_blank(),  # Remove panel background
                                plot.background = element_blank(),
                                panel.grid.minor = element_blank()
          )+
          theme(
            legend.position = "top",
            legend.direction = "horizontal"  # Arrange items in a row
          )+
          theme(text=element_text(size=21))
      
        p4<-ggplot(SpikeRatio, aes(x = Trial, y = SpikeRatio)) +
          geom_line() +
          labs(x = "Trial Number",
               y = "Modulation Index") +
          theme_bw()+
          ylim(-1,1)+
          geom_hline(yintercept = 0,color="blue",linetype="dashed")+
          scale_x_continuous(breaks = seq(1, length(trigger_times), by =2))+
          geom_hline(yintercept = mean(SpikeRatio$SpikeRatio),color="red",linetype="dashed")+
          coord_flip()+
          theme(                panel.background = element_blank(),  # Remove panel background
                                plot.background = element_blank())+
          theme(text=element_text(size=21))
        
        p5<-ggplot(k, aes(x = "", y = prop, fill = group)) +
          geom_bar(stat = "identity", width = 1, color = "white") +
          coord_polar("y", start = 0) +
          theme_void() + 
          geom_text(aes(y = ypos, label = paste0(round(prop), "%")), color = "white", size = 5,fontface="bold") +
          scale_fill_manual(values=c("red","black"))+
          theme(legend.title = element_blank())+
          theme(                panel.background = element_blank(),  # Remove panel background
                                plot.background = element_blank())+
          theme(
            legend.position = "top",
            legend.direction = "horizontal"  # Arrange items in a row
          )+
          theme(text=element_text(size=21))
        
       
        total_spike_counts_df$Bin_Start_Time<-as.numeric(total_spike_counts_df$Bin_Start_Time)
        data1<-total_spike_counts_df[total_spike_counts_df$Bin_Start_Time<(-0.5),]
        PreRate<-mean(data1$Total_Spike_Count)
        data2<-total_spike_counts_df[total_spike_counts_df$Bin_Start_Time>0,]
        
        
        ModulationIndex<-data.frame()
        for (i in 1:nrow(data2)){
          ModulationIndex[i,1]<-(data2[i,3]-PreRate)/(PreRate+data2[i,3])
        }
        
        ModulationIndex$StartTime<-data2$Bin_Start_Time
        colnames(ModulationIndex)<-c("MI","StartTime")
        
        p7<-ggplot(ModulationIndex,aes(x=StartTime,y=MI))+
          geom_line(alpha=0.5)+
          theme_bw()+
          geom_hline(yintercept = 0)+
          geom_smooth(method="lm",formula = y ~ poly(x, 15))+
          ylim(c(-1,1))+
          scale_x_continuous(breaks=seq(viewing_window[1],viewing_window[2],1))+
        xlab(c("Time since stimulus start (s)"))+
          ylab(c("Average modulation Index"))+
          theme(                panel.background = element_blank(),  # Remove panel background
                                plot.background = element_blank())+
          theme(text=element_text(size=21))
        
        p8<-ggplot(data = cluster_info_averaged_trials, aes(x = Rate_Before, y = Rate_After, color = Significant)) +
          geom_point() +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                             labels = c("TRUE" = "Modulated", "FALSE" = "Not Modulated")  # Change legend labels
          ) +
          theme_bw() +
          # geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
          xlab("Spontaneous FR (Hz)") + 
          ylab("Evoked FR (Hz)") +  
          labs(color=NULL)+
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
          theme(
            panel.background = element_blank(),  
            plot.background = element_blank())+
          theme(text=element_text(size=21))
        
         
  

            
        # p5.2<-ggplot(ModulationIndex,aes(x=StartTime,y=MI))+
        #   geom_line()+
        #   theme_bw()+
        #   geom_hline(yintercept = 0)+
        #   geom_smooth(method="lm",formula = y ~ poly(x, 15))+
        #   ylim(c(-1,1))+
        #   xlim(c(0,0.03))+
        #   xlab(c(""))+
        #   ylab(c(""))
        # 
        # 
        # 
        # 
        # p5<-p5.1+annotation_custom(ggplotGrob(p5.2),xmin=0,xmax=1.8,ymin=-1,ymax=-0.1)
        # 
        
        p6<-ggarrange(p5,p8,p2,p3,p1,p4,ncol=3, nrow=3,heights = c(1,1,1,1,1,1))
        
        list(p6 = p6,
             rasterplots=rasterplots,
             all_firing_rates_df = all_firing_rates_df, 
             total_spike_counts_df = total_spike_counts_df, 
             SpikeRatio = SpikeRatio)      })
      
      output$PSTH <- renderPlot({
        plot_data()$p6
      },height = 1400, width = 1800)
      
      output$Raster_Plots <- renderPlot({
        plot_data()$rasterplots
      },height = 1500, width = 1500)
      
      
      observeEvent(input$save_csv, {
        data <- plot_data()
        write.csv(data$all_firing_rates_df, "SpikeCount_FiringRate_PerBin_PerCluster.csv", row.names = FALSE)
        write.csv(data$total_spike_counts_df, "TotalSpikeCount_PerBin.csv", row.names = FALSE)
        write.csv(data$SpikeRatio, "ModulationIndex_PerTrial.csv", row.names = FALSE)
      })
      
      observeEvent(input$save_plot, {
        data <- plot_data()
        setwd("C:/Users/Antonio/Desktop")
        ggsave("PSTH.pdf", plot = data$p6, device = "pdf", width = 15, height = 15)
        ggsave("Rasters.pdf", plot = data$rasterplots, device = "pdf", width = 20, height = 15)
        
      })
      

  })
}
# Run the application
shinyApp(ui = ui, server = server)


