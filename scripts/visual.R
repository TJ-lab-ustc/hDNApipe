# 1. get user option 
args <- commandArgs(trailingOnly = TRUE)
visual_list <- c("category", "length", "patho", "circos", "go", "kegg", "ppi")
for (i in 1:length(visual_list)) {
    assign(visual_list[i], args[i])
}
num <- length(visual_list)
out_dir <- args[num+1]
SNP <- args[num+2]
INDEL <- args[num+3]
SV <- args[num+4]
CNV <- args[num+5]
bin <- as.numeric(args[num+6])
cadd_score <- as.numeric(args[num+7])

# 2. set variables
setwd(out_dir)
var_list <- c("SNP", "INDEL", "CNV", "SV")
chrs <- paste0('chr',c(1:22,'X','Y'))
for (var in var_list) {
    path <- get(var)
    if (path != "NA") {
        info=paste0("getting ", var, " vcf file......")
        print(info)
        vcf <- read.table(path)[, c(1,2,4,5,8)]
        colnames(vcf) <- c("chr", "site", "ref", "alt", "info")
        vcf <- vcf[vcf$chr %in% chrs,]
        assign(tolower(var), vcf)
    }
}
analysis=FALSE
if (go == "True" | kegg == "True" | ppi == "True")
    analysis=TRUE

# 3. import packages and define functions
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(stringr)
vep_annot_pre <- function(vcf_path, var_type) {
    info=paste0("getting ", toupper(var_type), " vcf annotation information......")
    print(info)

    title_line <- readLines(vcf_path)
    title_line <- title_line[grep("Format: ",title_line)]
    title_line <- strsplit(title_line, "Format: ")
    title_line <- strsplit(title_line[[1]][2], '"')
    title <- unlist(strsplit(title_line[[1]][1], split = "\\|"))

    vcf <- get(var_type)
    vcf_info <- separate(vcf, info, into = title, sep = "\\|", remove = TRUE)
    vcf_info <- vcf_info[, -c(1:4)]
    colnames(vcf_info) <- title

    return(vcf_info)
}

# 4. vcf pre-process
if (SV != 'NA') {
    print("getting SV vcf file info")
    x <- strsplit(sv$info, ";")
    sv$end_chr <- NA
    sv$end_pos <- NA
    sv$type <- NA
    for (i in 1:nrow(sv)) {
        type <- x[[i]][grep("^SVTYPE=", x[[i]])]
        sv[i,]$type <- sub("^SVTYPE=", "", type)
        if (sv[i,]$type == "BND") {
            string = sv[i,]$alt
            chr_index <- regexpr("chr", string)[1]
            colon_index <- regexpr(":", string)[1]
            chr <- substr(string, chr_index, colon_index - 1)
            pos <- str_extract(string, "(?<=:)\\d+")
            sv[i,]$end_chr <- chr
            sv[i,]$end_pos <- pos
            if (sv[i,]$end_chr == sv[i,]$chr)
            sv[i,]$type <- "INV"
        }
        else {
            end_pos <- x[[i]][grep("^END=", x[[i]])]
            sv[i,]$end_pos <- sub("^END=", "", end_pos) %>% as.integer()
            sv[i,]$end_chr <- sv[i,]$chr
        }
    }
    sv$end_pos <- as.numeric(sv$end_pos)
    sv$type <- factor(sv$type, levels = c("DUP", "DEL", "INS", "INV", "BND"))
}


# 5-1. category(raw)
# snp
if (category == "True" && SNP != 'NA') {

    # pre-process
    temp <- snp
    if (length(grep(",", snp$alt) != 0)) {
        print(paste(length(grep(",", snp$alt)), "sites with multiple SNP variants being ignored"))
        snp <- snp[grep(",", snp$alt, invert = T),]
    }

    # SNP classification
    snp$change_type1 <- NA
    snp$change_type2 <- 'Transversions'
    snp[(snp$ref == 'A' & snp$alt == 'T')|(snp$ref == 'T' & snp$alt == 'A'),]$change_type2 <- 'Transitions'
    snp[(snp$ref == 'C' & snp$alt == 'G')|(snp$ref == 'G' & snp$alt == 'C'),]$change_type2 <- 'Transitions'
    snp[(snp$ref == 'A' & snp$alt == 'T')|(snp$ref == 'T' & snp$alt == 'A'),]$change_type1 <- 'A>T|T>A'
    snp[(snp$ref == 'C' & snp$alt == 'G')|(snp$ref == 'G' & snp$alt == 'C'),]$change_type1 <- 'C>G|G>C'
    snp[(snp$ref == 'A' & snp$alt == 'C')|(snp$ref == 'T' & snp$alt == 'G'),]$change_type1 <- 'A>C|T>G'
    snp[(snp$ref == 'A' & snp$alt == 'G')|(snp$ref == 'T' & snp$alt == 'C'),]$change_type1 <- 'A>G|T>C'
    snp[(snp$ref == 'C' & snp$alt == 'A')|(snp$ref == 'G' & snp$alt == 'T'),]$change_type1 <- 'C>A|G>T'
    snp[(snp$ref == 'C' & snp$alt == 'T')|(snp$ref == 'G' & snp$alt == 'A'),]$change_type1 <- 'C>T|G>A'
    snp_type1 <- data.frame(table(snp$change_type1))
    colnames(snp_type1) <- c("type", "num")
    snp_type1$ratio <- round(100*snp_type1$num/sum(snp_type1$num), 2)
    snp_type2 <- data.frame(table(snp$change_type2))
    colnames(snp_type2) <- c("type", "num")
    snp_type2$ratio <- round(100*snp_type2$num/sum(snp_type2$num), 2)

    # plot
    print("plotting SNP base replacement......")
    p1 <- ggplot(data = snp_type1, aes(type, num)) +
        geom_bar(stat = 'identity', width = 0.6, fill='#3E549D') +
        labs(x = '', y = 'Count', title = "Distribution of single nucleotide mutation") +
        geom_text(aes(label = paste(ratio,"%")), vjust=-0.4, size = 3.5) +
        theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black", linewidth = 0.5),
                axis.title = element_text(size = 12), axis.text = element_text(size = 10),
                plot.title = element_text(family = "serif", size = 14, face = "bold", hjust = 0.4),
                legend.position = "none")
    ggsave('snp_BaseChange.pdf', plot=p1, width = 6, height = 4.5)
    print("*****SNP base replacement done*****")

    print("plotting SNP transitions/transversions......")
    ratio = snp_type2[snp_type2$type == "Transitions", ]$ratio/snp_type2[snp_type2$type == "Transversions", ]$ratio
    ratio = round(ratio, digits = 1)
    ratio_label = paste0("Ti/Tv ratio = ", ratio)
    p2 <- ggplot(data = snp_type2, aes(type, num)) +
        geom_bar(stat = 'identity', width = 0.6, fill = c("#DA3A2F", "#3E549D")) +
        geom_text(aes(label = paste(ratio,"%")), vjust=-0.4, size = 3.5) +
        labs(x = '', y = 'Count', title = ratio_label) +
        theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black", linewidth = 0.5),
                axis.title = element_text(size = 12), axis.text = element_text(size = 10),
                plot.title = element_text(family = "serif", size = 14, face = "bold", hjust = 0.4),
                legend.position = "none")
    ggsave('snp_Ti_Tv.pdf', plot=p2, width = 3, height = 4.5)
    print("*****SNP ti/tv done*****")
    
    snp <- temp
}

# indel
if (category =="True" && INDEL != 'NA') {
    # pre-process
    temp <- indel
    if (length(grep(",", indel$alt) != 0)) {
        print(paste(length(grep(",", indel$alt)),"sites with multiple INDEL variants being ignored"))
        indel <- indel[grep(",", indel$alt, invert = T),]
    }

    # count
    indel$len_change <- nchar(indel$alt) - nchar(indel$ref)
    indel$type <- ifelse(indel$len_change > 0, "insert", ifelse(indel$len_change < 0, "deletion", "else"))
    indel <- indel[indel$type != "else", ]
    indel$type <- factor(indel$type, levels = c("insert", "deletion"))

    # plot
    print("plotting INDEL type......")
    p1 <- ggplot(indel, aes(x = type, fill = type)) +
        geom_bar(stat = "count", width = 0.6) +
        scale_fill_manual(values=c("#DA3A2F", "#3E549D")) +
        labs(y="Count", x = "") +
        theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black", linewidth = 0.5),
                axis.title = element_text(size = 12), axis.text = element_text(size = 10),
                legend.position = "none")
    ggsave("indel_type.pdf", plot = p1, width = 3, height = 5)
    print("plotting INDEL type done*****")

    indel <- temp
}

# sv
if (category =="True" && SV != 'NA') {
    print("plotting SV type......")
    sv_types <- c("DUP", "DEL", "INS", "INV", "BND")
    sv_num <- c(sum(sv$type == "DUP"), sum(sv$type == "DEL"), sum(sv$type == "INS"), sum(sv$type == "INV"), sum(sv$type == "BND"))
    sv_info <- data.frame(sv_types,sv_num)
    sv_info$ratio <- paste(round(100*sv_info$sv_num/nrow(sv), 2), "%")
    pdf("sv_type.pdf", width = 6, height = 6)
    pie(x=sv_info$sv_num, labels = sv_info$ratio, col = c("#BAC4E2", "#8FA2D4", "#4472C4", "#315493", "#E6E6FA"), border = "white", radius = 0.9, main = "SV type", family = "Times")
    legend("topright", c("DUP", "DEL", "INS", "INV", "BND"), cex = 1, fill = c("#BAC4E2", "#8FA2D4", "#4472C4", "#315493", "#E6E6FA"), box.lty=0 )
    dev.off()
    print("*****SV type done*****")
}

# 5-2. length(raw)
# indel
if (length =="True" && INDEL != 'NA') {
    print("plotting INDEL length......")
    temp <- indel
    indel$len_change <- nchar(indel$alt) - nchar(indel$ref)
    indel$type <- ifelse(indel$len_change > 0, "insert", ifelse(indel$len_change < 0, "deletion", "else"))
    indel <- indel[indel$type != "else", ]
    indel$type <- factor(indel$type, levels = c("insert", "deletion"))
    p2 <- ggplot(indel, aes(x = type, y = abs(len_change), color = type)) +
        geom_jitter(size = 0.5) +
        scale_color_manual(values=c("#DA3A2F", "#3E549D")) +
        labs(x = "", y = "Length change") +
        theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black", linewidth = 0.5),
                axis.title = element_text(size = 12), axis.text = element_text(size = 10),
                legend.position = "none")
    ggsave("indel_length.pdf", plot = p2, width = 3, height = 5)
    indel <- temp
    print("*****INDEL length done*****")
}

# sv
if (length =="True" && SV != 'NA') {
    print("plotting SV length......")
    temp <- sv
    sv <- sv[sv$type != "BND", ]
    sv$len_log <- NA
    sv$len_log <- log10(abs(sv$end_pos - sv$site))
    for (i in 1:nrow(sv)) {
        if (sv$type[i] == "INS")
            sv$len_log[i] <- log10(nchar(sv$alt[i]) - nchar(sv$ref[1]))
    }
    sv <- sv[sv$len_log != -Inf, ]
    p1 <- ggplot(sv, aes(x=len_log,fill=type)) +
        geom_histogram(binwidth = 0.5,position="stack", color = "white", linewidth = 0.1, center = 0.25, alpha = 0.9) +
        scale_fill_manual(values=c("#BAC4E2", "#8FA2D4", "#4472C4", "#315493")) +
        theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black",linewidth = 0.5),
                axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10),
                plot.title = element_text(family = "serif", size = 14, face = "bold", hjust = 0.5),
                legend.text = element_text(size = 10), legend.title = element_text(size = 10, face = "bold"))+
        labs(x = "Length(log10)", y="Number of SV", title = "SV length distribution", fill = "Type") +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), limits = c(1, ceiling(max(sv$len_log))), breaks = seq(0,ceiling(max(sv$len_log)), 1))
    ggsave("sv_length.pdf", plot = p1, width = 6, height = 4.5)
    sv <- temp
    print("*****SV length done*****")
}

# 5-3. Circos (raw)
if (circos == 'True' && SNP != 'NA' && INDEL != 'NA' && SV != 'NA') {
    
    print("plotting Circos plot......")
    
    # packages
    library(circlize)
    library(ComplexHeatmap)
    library(gridBase)

    # snp
    snp$site <- as.numeric(snp$site)
    snp$chr <- as.factor(snp$chr)

    # chromosome information
    chr_path <- paste0(out_dir, "/chr_len.txt")
    chrlen <- read.table(chr_path,header = F)
    chrlen$end <- chrlen$V2
    chrlen$V2 <- 1
    
    # snp
    for (chr in chrs) {
        subset <- snp[snp$chr == chr, ]
        start_pos <- 1
        end_pos <- bin
        region_num <- floor(chrlen[chrlen$V1 == chr, ]$end/bin)
        snp_info <- data.frame(chr=chr, start=0, end=0, snp_num=0)[-1,]
        for (i in 1:region_num) {
            snp_info[i,]$start <- start_pos + (i-1)*bin
            snp_info[i,]$end <- end_pos + (i-1)*bin
            snp_info[i,]$snp_num <- length(which(subset$site <= snp_info[i,]$end & subset$site >= snp_info[i,]$start))
        }
        if (chrlen[chrlen$V1 == chr,]$end %% bin != 0) {
            snp_info[(i+1),]$start <- snp_info[i,]$end + 1
            snp_info[(i+1),]$end <- chrlen[chrlen$V1 == chr,]$end
            snp_info[(i+1),]$snp_num <- length(which(subset$site > snp_info[i,]$end))
        }
        snp_info$chr <- chr
        assign(paste0('snp_',chr), snp_info)
    }
    snp_info <- data.frame(chr=0, start=0, end=0, snp_num=0)[-1,]
    for (chr in chrs) {
        x <- get(paste0('snp_',chr))
        snp_info <- rbind(snp_info,x)
    }
    snp_legend <- Legend(at = c("SNP"), type = "lines", background = "white", size = 2,
                        legend_gp = gpar(col = "#8C9299"), title_position = "topleft", 
                        title = "SNP" )

    # indel
    indel$len_change <- nchar(indel$alt) - nchar(indel$ref)
    indel$type <- 'else'
    indel[indel$len_change>0,]$type <- 'insert'
    indel[indel$len_change<0,]$type <- 'deletion'
    for (chr in chrs) {
        subset <- indel[indel$chr==chr, ]
        subset <- subset[grep(",", subset$ref, invert = T), ]
        start_pos <- 1
        end_pos <- bin
        region_num <- floor(chrlen[chrlen$V1==chr,]$end/bin)
        indel_info <- data.frame(chr=chr, start=0, end=0, insert_num=0, deletion_num=0)[-1,]
        for (i in 1:region_num){
            indel_info[i,]$start <- start_pos+(i-1)*bin
            indel_info[i,]$end <- end_pos+(i-1)*bin
            indel_info[i,]$insert_num <- length(which(subset$site<=indel_info[i,]$end & subset$site>=indel_info[i,]$start & subset$type=='insert'))
            indel_info[i,]$deletion_num <- length(which(subset$site<=indel_info[i,]$end & subset$site>=indel_info[i,]$start & subset$type=='deletion'))
        }
        if (chrlen[chrlen$V1==chr,]$end %% bin != 0){
            indel_info[(i+1),]$start <- indel_info[i,]$end + 1
            indel_info[(i+1),]$end <- chrlen[chrlen$V1==chr,]$end
            indel_info[(i+1),]$insert_num <- length(which(subset$site>indel_info[i,]$end & subset$type=='insert'))
            indel_info[(i+1),]$deletion_num <- length(which(subset$site>indel_info[i,]$end & subset$type=='deletion'))
        }
        indel_info$chr <- chr
        assign(paste0('indel_', chr), indel_info)
    }
    indel_info <- data.frame(chr=0, start=0, end=0, indel_num=0)[-1,]
    for (chr in chrs) {
        x <- get(paste0('indel_', chr))
        indel_info <- rbind(indel_info, x)
    }
    indel_info$deletion_num <- -indel_info$deletion_num
    indel_legend = Legend(at = c("Insert", "Deletion"), type = "points", pch = NA, background = c("#EE6842","#6BB3DC"),
                            title_position = "topleft", row_gap = unit(1, "mm"),
                            title = "INDEL")

    # clean
    x <- ls()
    variables_to_remove <- grep("^snp_chr|^indel_chr", x, value = TRUE)
    rm(list = variables_to_remove)

    # cnv
    if (CNV != 'NA' ) {
        x <- strsplit(cnv$info, ";")
        cnv$end_pos <- NA
        cnv$fc_log <- NA
        for (i in 1:nrow(cnv)) {
            end_pos <- x[[i]][grep("END=", x[[i]])]
            cnv[i,]$end_pos <- sub("^END=", "", end_pos) %>% as.integer()
            fc_log <- x[[i]][grep("FOLD_CHANGE_LOG=", x[[i]])]
            cnv[i,]$fc_log <- sub("^FOLD_CHANGE_LOG=", "", fc_log) %>% as.numeric()
        }
        cnvmin <- min(cnv$fc_log)
        cnvmax <- max(cnv$fc_log)
        cnv_legend <- Legend(at = c("Duplication","Deletion"), type = "points", pch = NA, background = c("#ed1c21","#6ACDE0"),
                            title_position = "topleft", row_gap = unit(1, "mm"),
                            title = "CNV")
    }
    

    # sv
    sv_dup <- sv[sv$type == "DUP", c("chr", "site", "end_pos", "type")]
    sv_del <- sv[sv$type == "DEL", c("chr", "site", "end_pos", "type")]
    sv_inv <- sv[sv$type == "INV",][, c("chr", "site", "end_chr", "end_pos", "type")]
    sv_bnd <- sv[sv$type == "BND",][, c("chr", "site", "end_chr", "end_pos", "type")]
    sv_legend <- Legend(at = c("Duplication", "Deletion", "Inversion", "Translocation"),type = c("points", "points", "lines", "lines"),
                        pch = NA, background = c("#cf3a40", "#B2DA70", "white", "white"), legend_gp = gpar(col=c(NA, NA, "#8C9299", "#E9BF44")), row_gap = unit(1, "mm"),
                        title_position = "topleft", title = "SV")

    # legend
    if (CNV != 'NA')
        lgd_list_vertical <- packLegend(snp_legend, indel_legend, cnv_legend, sv_legend, direction = "vertical")
    else
        lgd_list_vertical <- packLegend(snp_legend, indel_legend, sv_legend, direction = "vertical")

    # plot
    pdf("circos.pdf", width = 10, height = 8)
    plot.new()
    circle_size = unit(1, "snpc")
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circos.par(start.degree = 90)
    circos.initializeWithIdeogram(species = 'hg38', track.height = 0.05)
    # snp
    circos.genomicTrackPlotRegion(snp_info, bg.border="#BFBFBF", track.height = 0.1, panel.fun = function(region, value, ...){
        circos.genomicLines(region, value, ytop.column = 1, col="#8C9299", border = "#8C9299", lwd = 2, ...)
    })
    # indel
    circos.genomicTrackPlotRegion(indel_info, bg.border="#BFBFBF", track.height=0.2, panel.fun = function(region, value, ...){
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,col = "#EE6842", border = "#EE6842", ...)
        circos.genomicRect(region, value, ytop.column = 2, ybottom = 0,col = "#6BB3DC", border = "#6BB3DC", ...)
    })
    # cnv
    if (CNV != 'NA') {
        circos.genomicTrackPlotRegion(cnv[,c(1,2,6,7)], bg.border="#BFBFBF", ylim=c(cnvmin,cnvmax), track.height=0.1, panel.fun = function(region, value, ...){
        col_val = ifelse(value[[1]] > 0, "#ed1c21", "#6ACDE0")
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,col = col_val, border = col_val, ...)
        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#DDDDDD")})
    }

    # sv
    set_track_gap(mm_h(0))
    circos.genomicTrackPlotRegion(sv_dup, bg.border= NA, ylim=c(0,1), track.height=0.05, panel.fun = function(region, value, ...){
        circos.genomicRect(region, value, col= "#cf3a40", border = "#cf3a40", ...)
    })
    circos.genomicTrackPlotRegion(sv_del, bg.border= NA, ylim=c(0,1), track.height=0.05, panel.fun = function(region, value, ...){
        circos.genomicRect(region, value, col= "#B2DA70", border = "#B2DA70", ...)
    })
    circos.genomicLink(sv_inv[,c(1,2,2)], sv_inv[,c(3,4,4)], col = "#8C9299", lwd = 0.5)
    circos.genomicLink(sv_bnd[,c(1,2,2)], sv_bnd[,c(3,4,4)], col = "#E9BF44", lwd = 0.5)
    # legend
    upViewport()
    draw(lgd_list_vertical, x = circle_size,  just = c("left"))
    dev.off()

    print("*****Circos plot done*****")
}

# 5-4. pathogenicity and analysis of short variants (annotated)
# snp (SNFT/PolyPhen2 annotated)
genelist <- NULL
if ((patho =="True" | analysis) & SNP != 'NA') {
    snp_info <- vep_annot_pre(SNP, "snp")
    if("SYMBOL" %in% colnames(snp_info)) {
        pathogenicity <- snp_info[(snp_info$SIFT != "" & snp_info$PolyPhen != ""), c("Gene","SYMBOL", "SIFT", "PolyPhen")]
    } else {
        pathogenicity <- snp_info[(snp_info$SIFT != "" & snp_info$PolyPhen != ""), c("Gene","SIFT", "PolyPhen")]
    }
    pathogenicity$SIFT <- gsub(".*\\(","",pathogenicity$SIFT)
    pathogenicity$SIFT <- gsub(")","",pathogenicity$SIFT) %>% as.numeric()
    pathogenicity$PolyPhen <- gsub(".*\\(","",pathogenicity$PolyPhen)
    pathogenicity$PolyPhen <- gsub(")","",pathogenicity$PolyPhen) %>% as.numeric()
    pathogenicity$sig <- ifelse (
        pathogenicity$SIFT < 0.05 & pathogenicity$PolyPhen > 0.908, "Deleterious",
        ifelse (pathogenicity$SIFT >= 0.05 & pathogenicity$PolyPhen <= 0.446, "Benign", "Possible")
        )
    pathogenicity$sig <- factor(pathogenicity$sig, levels = c("Benign", "Possible", "Deleterious"))
    write.csv(pathogenicity, "snp_pathogenicity.csv", row.names = FALSE)
    # plot pathogenicity
    if (patho == "True") { 
        print("plotting SNP SIFT/PolyPhen2......")
        p1 <- ggplot(pathogenicity, aes(x = PolyPhen, y = SIFT, color = sig)) +
            geom_point(size = 0.5) +
            geom_vline(xintercept = 0.446, linetype = "dashed", color = "#cde2d0", linewidth = 0.7) +
            geom_vline(xintercept = 0.908, linetype = "dashed", color = "#cde2d0", linewidth = 0.7) +
            geom_hline(yintercept = 0.05, linetype = "dashed", color = "#cde2d0", linewidth = 0.7) +
            labs(title = "Pathogenicity prediction") +
            scale_color_manual(values = c("grey", "#85c8ed", "#ff8b73")) +
            theme(panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.5),
                axis.text = element_text(family = "serif", size = 10),
                axis.title = element_text(family = "serif", face = "bold", size = 12),
                plot.title = element_text(family = "serif", face = "bold", hjust = 0.5, size = 14),
                legend.title = element_text(family = "serif", face = "bold", size = 12),
                legend.text = element_text(family = "serif", size =10),
                legend.key = element_blank()) +
            guides(color=guide_legend(override.aes = list(size=1,alpha=1)))
        ggsave("snp_pathogenicity.pdf", plot = p1, width = 6, height = 5)
        print("*****SNP SIFT/PolyPhen2 done*****")
    }
    if (analysis)
        genelist<-c(genelist, pathogenicity[pathogenicity$sig!="Benign", "Gene"])
}

# indel (CADD annotated)
if ((patho =="True" | analysis) & INDEL != 'NA') {
    indel_info <- vep_annot_pre(INDEL, "indel")
    if ("SYMBOL" %in% colnames(indel_info)) {
      cadd_table <- indel_info[indel_info$CADD_PHRED != "" | indel_info$CADD_RAW != "", c("Gene", "SYMBOL", "CADD_PHRED", "CADD_RAW")]
    } else {
      cadd_table <- indel_info[indel_info$CADD_PHRED != "" | indel_info$CADD_RAW != "", c("Gene","CADD_PHRED", "CADD_RAW")]
    }
    cadd_table$CADD_PHRED <- as.numeric(cadd_table$CADD_PHRED)
    cadd_table$CADD_RAW <- as.numeric(cadd_table$CADD_RAW)
    write.csv(cadd_table, "indel_cadd.csv", row.names = FALSE)
    if (patho == "True") {  
        print("plotting INDEL CADD......")
        cadd_phred <- indel_info[indel_info$CADD_PHRED != "" , "CADD_PHRED"] %>% as.data.frame()
        cadd_phred$group <- "PHRED"
        cadd_raw <- indel_info[indel_info$CADD_RAW != "" , "CADD_RAW"] %>% as.data.frame()
        cadd_raw$group <- "RAW"
        cadd <- rbind(cadd_phred, cadd_raw)
        colnames(cadd)[1] <- "score"
        cadd$score <- as.numeric(cadd$score)
        p1 <- ggplot(data = cadd, aes(x = group, y = score, fill = group)) +
            geom_violin(color = "white") +
            geom_boxplot(width = 0.05, outlier.size = 0.5) +
            scale_fill_manual(values = c("#55B5E8", "#E59F01")) +
            theme(panel.grid=element_blank(), panel.background = element_rect(fill = 'white'), axis.line = element_line(color = "black", linewidth = 0.5),
                    axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10),
                    legend.position = "none")
            ggsave("indel_cadd.pdf", plot = p1, width = 4, height = 4)
        print("*****INDEL CADD done*****")
    }
    if (analysis) 
        genelist <- c(genelist, cadd_table[cadd_table$CADD_PHRED > cadd_score, "Gene"])
}


# 6. analysis
if (analysis) {
    
    genelist <- genelist[!duplicated(genelist)]
    write.csv(genelist, "genelist.csv", row.names = FALSE)

    library(clusterProfiler)
    library(org.Hs.eg.db)
    gene_entre <- bitr(genelist, fromType = "ENSEMBL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

    if (go == "True") {
        ego <- enrichGO(gene = gene_entre$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "all",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
        ego_res <- ego@result
        ego_res <- ego_res[order(ego_res$qvalue,decreasing=F), ]
        write.csv(ego_res, "GO.csv")
        ego_res <- ego_res[ego_res$pvalue < 0.05, ]
        if (nrow(ego_res) == 0) {
            writeLines("Warning: No items were enriched in the GO analysis with a p-value under 0.05.\nSkipping plotting.")
        } else {
            for (i in 1:nrow(ego_res)) {
                ego_res$Ratio[i] <- as.numeric(strsplit(ego_res$GeneRatio[i], split = "/")[[1]][1])/as.numeric(strsplit(ego_res$BgRatio[i], split = "/")[[1]][1])
            }
            ego_res$Category <- ego_res$ONTOLOGY
            ego_res$Category <- gsub("CC", "cellular component", ego_res$Category)
            ego_res$Category <- gsub("BP", "biological process", ego_res$Category)
            ego_res$Category <- gsub("MF", "molecular function", ego_res$Category)
            if (nrow(ego_res) > 30) ego_res <- ego_res[1:30, ]
            ego_res <- ego_res[order(ego_res$qvalue,decreasing=T), ]
            level <- ego_res$Description
            ego_res$Description <- factor(ego_res$Description, levels = level)
            p1 <- ggplot(ego_res, aes(x = Ratio,y = Description, shape = Category)) +
                geom_point(aes(size = Count, color = qvalue)) +
                theme_bw() +
                labs(y = "", x = "RichFactor", title = "GO enrichment") +
                scale_color_gradient(low = "red", high = "blue") +
                theme(axis.text = element_text(face = "bold", size = 10),
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(face = "bold", size = 12),
                    plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
            ggsave("GO.pdf", plot = p1, width = 8, height = (3 + 0.2*nrow(ego_res)))
            print("*****GO Analysis done*****")
        }
    }


    if (kegg == "True") {
        kk <- enrichKEGG(gene = gene_entre$ENTREZID, organism ="hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        kk_res <- kk@result
        kk_res <- kk_res[order(kk_res$qvalue,decreasing=F), ]
        write.csv(kk_res, "KEGG.csv")
        kk_res <- kk_res[kk_res$pvalue < 0.05, ]
        if (nrow(kk_res) == 0) {
            writeLines("Warning: No items were enriched in the KEGG analysis with a p-value under 0.05.\nSkipping plotting.")
        } else {
            for (i in 1:nrow(kk_res)) {
                kk_res$Ratio[i] <- as.numeric(strsplit(kk_res$GeneRatio[i],split = "/")[[1]][1])/as.numeric(strsplit(kk_res$BgRatio[i],split = "/")[[1]][1])
            }
            if (nrow(kk_res) > 30) kk_res <- kk_res[1:30, ]
            kk_res <- kk_res[order(kk_res$qvalue,decreasing=T), ]
            level <- kk_res$Description
            kk_res$Description <- factor(kk_res$Description, levels = level)
            p1 <- ggplot(kk_res, aes(x=Ratio,y=Description)) +
                geom_point(aes(size=Count,color=qvalue)) +
                theme_bw() +
                labs(y="", x="RichFactor", title = "KEGG enrichment") +
                scale_color_gradient(low="red",high="blue") +
                theme(axis.text = element_text(face = "bold", size = 10),
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(face = "bold", size = 12),
                    plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
            ggsave("KEGG.pdf", plot = p1, width = 8, height = (4+0.2*nrow(kk_res)))
            print("*****KEGG Analysis done*****")
    }
    }

    if (ppi == "True") {
        print("running PPI Analysis......")
        library(STRINGdb)
        library(ggraph)
        genelist1 <- bitr(gene_entre$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
        string_db <- STRINGdb$new(version="11.5", species=9606, network_type="physical")
        data_mapped <- string_db$map(genelist1, my_data_frame_id_col_names = "ENTREZID", removeUnmappedRows = TRUE)
        pdf("PPI_raw.pdf", width = 10, height = 10)
        string_db$plot_network(data_mapped$STRING_id)
        dev.off()
        data_links <- data_mapped$STRING_id %>% string_db$get_interactions() %>% unique() 
        data_links$from <- data_mapped[match(data_links$from, data_mapped$STRING_id), "SYMBOL"]
        data_links$to <- data_mapped[match(data_links$to, data_mapped$STRING_id), "SYMBOL"]
        pdf("PPI_processed.pdf", width = 6, height = 5.5)
        data_links %>% tidygraph::as_tbl_graph() %>%
            ggraph(layout = "graphopt")+
            geom_edge_fan(color = "grey")+
            geom_node_point(size=5, color="lightblue", alpha = 1)+
            geom_node_text(aes(label=name), repel = T)+
            theme_void() +
            labs(title = "Protein-protein interaction") +
            theme(plot.title = element_text(family = "serif", face = "bold", size = 12, hjust = 0.5))
        dev.off()
        print("*****PPI Analysis done*****")
    }

}

writeLines("*****Plotting done.*****\n\n")
