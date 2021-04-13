#Jenny Smith
# March 30, 2021
#Purpose: Create donut plots for SNVs types of interest.

#' Create a donut plot
#'
#' @param variants.df is the filtered and re-formatted dataframe from filtered VCF from filter_vcf() function.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' vcf <- vcfR::read.vcfR("/path/to/vcf")
#' vcf.df <- cbind(as.data.frame(vcfR::getFIX(vcf)), vcfR::INFO2df(vcf))
#' vcf.s4 <- VariantAnnotation::readVcf("/path/to/vcf")
#' variants.df <- filter_vcf(vcf.df,vcf.s4)
#'
#' donut_plot(filtered.vcf=variants.df)
#' }
#'
#' @import dplyr
donut_plot <- function(filtered.vcf){

  #Calculate number of variant types
  donut.df <- filtered.vcf %>%
    group_by(Annotation_1) %>%
    summarize(N=n()) %>%
    mutate_at(vars(Annotation_1), ~case_when(
      grepl("5_prime_UTR", .) ~ "5_prime_UTR_variant",
      grepl("missense_variant", .) ~ "missense_variant",
      grepl("splice", .) ~ "splice_site",
      TRUE ~ .)) %>%
    group_by(Annotation_1) %>%
    summarise(N=sum(N)) %>%
    ungroup() %>%
    mutate(Percent=round(N/sum(N) *100, digits = 2)) %>%
    arrange(N) %>%
    mutate(lab.pos = cumsum(Percent)-.5*Percent) %>%
    mutate(Annotation_1 = factor(Annotation_1, levels=rev(unique(Annotation_1))))

  #subset dataset if necessary
  GT.5pct <- filter(donut.df, Percent > 5)
  LE.5pct <- filter(donut.df, Percent <= 5)

  #create the donut plot with ggplot2
  donut.types <- ggplot2::ggplot(data = donut.df,
                        aes(x = 2, y = Percent, fill = Annotation_1))+
    ggplot2::geom_col(width = 1.0, color="black", size=0.25) +
    # scale_x_discrete(limits = c(" ", 2)) +
    ggplot2::xlim(0.3, 2.5) +
    ggplot2::annotate(geom="label",
             x=2.0, y=GT.5pct$lab.pos,
             label=paste0(GT.5pct$Percent,"%"),
             size=3) +
    ggplot2::labs(title="Variant Locations within Genes") +
    ggplot2::coord_polar("y", start=1) +
    ggplot2::scale_fill_brewer(palette="Paired") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "left",
          legend.title = element_blank(),
          legend.text = element_text(size=8))

  #if there are some categories with low frequency, add thier labels using ggrepel
  if(nrow(LE.5pct) > 0){
    donut.types <- donut.types +
      ggrepel::geom_text_repel(data=LE.5pct,
                               aes(x = 2.5, y = lab.pos,
                                   label=paste0(Percent,"%")),
                               nudge_x = 0.25,
                               segment.size = .5,
                               min.segment.length = 0.1,
                               size=3,
                               show.legend = FALSE,
                               inherit.aes = F)
  }

  return(donut.types)
}
