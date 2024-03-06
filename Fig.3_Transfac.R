library(sqldf)

setwd('G:/Мой диск/Bioinform/Dissertation/Доп.мат/Plots/tfbs')

COL=pal_locuszoom()(4)



svg("G:/Мой диск/Bioinform/Dissertation/Графики/20.TRANSFAC_motifs.svg", width = 10, height = 7)
par(mfrow=c(2,2))

                    ## AP2_EREBP-related_factors
FILES=list.files(pattern = 'AP2_EREBP-related')
plot(0, xlim=c(-600,200), ylim=c(0.001,0.009), 
     pch = '',
     ylab ='Частота' , xlab = 'Расстояние от сайта начала транскрипции, п.н.',
     main='AP2/EREBP-related factors')

for(i in 1:length(FILES)){
  f=FILES[i]
  Data=read.table(f, sep='\t')
  names(Data)=c('Gene','Matrix','START','SEQ','END','CoreSim','MatSim')
  # names(Data)=c('Gene','Start','End')  
od1=sqldf('select distinct Gene, round(((START+END)/2-400)/5,0)*5 from Data')
names(od1)[2]='POS'
  #od1$POS=round(0.5*(od1$START+od1$END)/20,0)*20
  od2=aggregate(od1$Gene, by=list(od1$POS), FUN='length')
  S=sum(od2$x)
  lines(od2$Group.1, od2$x/S, col=COL[i], lwd=2)
}

abline(v=0)

labs=c(paste('L. sibirica, ', nrow(read.table(FILES[1], sep='\t')), 'occurances'),
       paste('P. abies, ', nrow(read.table(FILES[2], sep='\t')), 'occurances'),
       paste('P. glauca, ', nrow(read.table(FILES[3], sep='\t')), 'occurances'),
       paste('P. taeda, ', nrow(read.table(FILES[4], sep='\t')), 'occurances'))
legend('topleft', legend=labs, lwd=2, col=COL, text.font=1, bty='n', cex=0.9, y.intersp = 0.8)


                ## Homeo_domain
FILES=list.files(pattern = 'Homeo_domain')
  plot(0, xlim=c(-600,200), ylim=c(0.001,0.009), 
     pch = '',
     ylab ='Частота' , xlab = 'Расстояние от сайта начала транскрипции, п.н.',
     main='Homeo domain')

  for(i in 1:length(FILES)){
  f=FILES[i]
  Data=read.table(f, sep='\t')
  names(Data)=c('Gene','Matrix','START','SEQ','END','CoreSim','MatSim')
  # names(Data)=c('Gene','Start','End')  
od1=sqldf('select distinct Gene, round(((START+END)/2-400)/5,0)*5 from Data')
names(od1)[2]='POS'
  #od1$POS=round(0.5*(od1$START+od1$END)/20,0)*20
  od2=aggregate(od1$Gene, by=list(od1$POS), FUN='length')
  S=sum(od2$x)
  lines(od2$Group.1, od2$x/S, col=COL[i], lwd=2)
  }
  
  abline(v=0)
  
labs=c(paste('L. sibirica, ', nrow(read.table(FILES[1], sep='\t')), 'occurances'),
       paste('P. abies, ', nrow(read.table(FILES[2], sep='\t')), 'occurances'),
       paste('P. glauca, ', nrow(read.table(FILES[3], sep='\t')), 'occurances'),
       paste('P. taeda, ', nrow(read.table(FILES[4], sep='\t')), 'occurances'))
legend('topleft', legend=labs, lwd=2, col=COL, text.font=1, bty='n', cex=1, y.intersp = 0.8)
  
              ## Heat Schock factors
FILES=list.files(pattern = 'HSF')
  plot(0, xlim=c(-600,200), ylim=c(0.001,0.009), 
     pch = '',
     ylab ='Частота' , xlab = 'Расстояние от сайта начала транскрипции, п.н.',
     main='Heat shock factors')

  for(i in 1:length(FILES)){
  f=FILES[i]
  Data=read.table(f, sep='\t')
  names(Data)=c('Gene','Matrix','START','SEQ','END','CoreSim','MatSim')
  # names(Data)=c('Gene','Start','End')  
od1=sqldf('select distinct Gene, round(((START+END)/2-400)/5,0)*5 from Data')
names(od1)[2]='POS'
  #od1$POS=round(0.5*(od1$START+od1$END)/20,0)*20
  od2=aggregate(od1$Gene, by=list(od1$POS), FUN='length')
  S=sum(od2$x)
  lines(od2$Group.1, od2$x/S, col=COL[i], lwd=2)
  }
  
  abline(v=0)

  labs=c(paste('L. sibirica, ', nrow(read.table(FILES[1], sep='\t')), 'occurances'),
       paste('P. abies, ', nrow(read.table(FILES[2], sep='\t')), 'occurances'),
       paste('P. glauca, ', nrow(read.table(FILES[3], sep='\t')), 'occurances'),
       paste('P. taeda, ', nrow(read.table(FILES[4], sep='\t')), 'occurances'))
legend('topleft', legend=labs, lwd=2, col=COL, text.font=1, bty='n', cex=1, y.intersp = 0.8)
  
          ## MYB factors
FILES=list.files(pattern = 'MYB')
  plot(0, xlim=c(-600,200), ylim=c(0.001,0.009), 
     pch = '',
     ylab ='Частота' , xlab = 'Расстояние от сайта начала транскрипции, п.н.',
     main='MYB')

  for(i in 1:length(FILES)){
  f=FILES[i]
  Data=read.table(f, sep='\t')
  names(Data)=c('Gene','Matrix','START','SEQ','END','CoreSim','MatSim')
  # names(Data)=c('Gene','Start','End')  
od1=sqldf('select distinct Gene, round(((START+END)/2-400)/5,0)*5 from Data')
names(od1)[2]='POS'
  #od1$POS=round(0.5*(od1$START+od1$END)/20,0)*20
  od2=aggregate(od1$Gene, by=list(od1$POS), FUN='length')
  S=sum(od2$x)
  lines(od2$Group.1, od2$x/S, col=COL[i], lwd=2)
  }
  
  abline(v=0)

  labs=c(paste('L. sibirica, ', nrow(read.table(FILES[1], sep='\t')), 'occurances'),
       paste('P. abies, ', nrow(read.table(FILES[2], sep='\t')), 'occurances'),
       paste('P. glauca, ', nrow(read.table(FILES[3], sep='\t')), 'occurances'),
       paste('P. taeda, ', nrow(read.table(FILES[4], sep='\t')), 'occurances'))
legend('topleft', legend=labs, lwd=2, col=COL, text.font=1, bty='n', cex=1, y.intersp = 0.8)

dev.off()
