if (file.exists('~/.Rprofile')) {
  source('~/.Rprofile')
}

library(knitr)
library(magrittr)
options(knitr.table.format = 'markdown')

clean_post <- function(prefix, fig_path) {
  paste0("^", prefix) %>% 
    list.files(path = fig_path, pattern = ., full.names = TRUE) %>% 
    unlink
}

post_title <- function(source_rmd) {
  basename(source_rmd) %>% 
  sub("^[0-9]+-[0-9]+-[0-9]+-", "", .) %>% 
    sub("\\.rmd$", "", ., ignore.case = TRUE)
}

is_published <- function(source_rmd) {
  published <- readLines(source_rmd, warn = FALSE, n = 25) %>% 
    tolower %>% 
    grep("published", ., value = TRUE) %>% 
    sapply(function(x) grepl('true', x, fixed = TRUE))
  if (length(published) > 0 && !published[[1]]) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

knit_fig_path <- function(input, output, fig_path, cache_path, ...) {
  knitr::opts_chunk$set(fig.path = fig_path, cache.path = cache_path)
  knitr::knit(input = input, output = output, ...)
  invisible()
}

hook_plot_md_link <- function(x, options) {
  if ('img.link' %in% names(options) && options$img.link) {
    paste0(
      sprintf('<a href="%s">', paste0(opts_knit$get('base.url'), x)),
      knitr::hook_plot_md(x, options), 
      '</a>')
  } else {
    knitr::hook_plot_md(x, options)
  }
}

# From: http://chepec.se/2014/07/16/knitr-jekyll.html
knit_post <- function(source_rmd = '', img.dev = "svg", overwrite = FALSE, clean = FALSE) {
  # local directory of jekyll site
  site_path <- '/Users/mes335/projects/mstrimas.github.io/'
  # rmd directory (relative to base)
  rmd_path <- file.path(site_path, '_source/')
  # relative directory of figures
  fig_url <- 'figures/'
  # directory for markdown output
  posts_path <- file.path(site_path, '_posts/')
  
  # knitr options
  #knitr::render_jekyll(highlight = 'pygments')
  knitr::render_markdown(strict = F)
  knitr::knit_hooks$set(plot = hook_plot_md_link)
  knitr::opts_knit$set(
    base.url='/',
    base.dir=site_path)
  knitr::opts_chunk$set(
    fig.path=fig_url,
    fig.align='center',
    dpi=96,
    fig.width=480/96,
    fig.height=480/96,
    dev=img.dev,
    comment='#>',
    tidy=FALSE,
    cache=FALSE,
    error=TRUE,
    warning=FALSE,
    message=FALSE,
    collapse=TRUE)
  
  if (source_rmd == '') {
    files_rmd <- dplyr::data_frame(rmd = list.files(
      path = rmd_path,
      full.names = TRUE,
      pattern = '\\.rmd$',
      ignore.case = TRUE,
      recursive = FALSE))
    # create list of files to knit
    file_df <- files_rmd %>% 
      #dplyr::group_by(rmd) %>% 
      dplyr::mutate(
        base_name = gsub('\\.rmd$', '', basename(rmd)),
        title = post_title(rmd),
        md = file.path(posts_path, paste0(base_name, '.md')),
        fig_path = file.path(fig_url, paste0(title, '_')),
        cache_path = file.path(rmd_path, "cache", paste0(title, '/')),
        md_exists = file.exists(md),
        published = sapply(rmd, is_published),
        md_render = ifelse(published & (overwrite | !md_exists), TRUE, FALSE)) %>%
      dplyr::filter(md_render)
    
    if (nrow(file_df) == 0) {
      return(invisible())
    }
    
    # clean
    if (clean) {
      file_df %>% 
        dplyr::select(rmd, title) %>% 
        plyr::a_ply(1, transform, clean_post(
          prefix = title, fig_path = file.path(site_path, fig_url)))
    }
    
    # knit
    file_df %>% 
      dplyr::select(rmd, md, fig_path, cache_path) %>% 
      plyr::a_ply(1, transform, 
                  knit_fig_path(
                    input = rmd, output = md, fig_path = fig_path, 
                    cache_path = cache_path, quiet = TRUE))
  } else {
    source_rmd <- file.path(rmd_path, source_rmd)
    stopifnot(file.exists(source_rmd))
    base_name <- gsub('\\.rmd$', '', basename(source_rmd))
    title <- post_title(source_rmd)
    md_file <- file.path(posts_path, paste0(base_name, '.md'))
    fig_path <- file.path(fig_url, paste0(title, '_'))
    cache_path <- file.path(rmd_path, 'cache', paste0(title, '/'))
    if (clean) {
      clean_post(prefix = title, fig_path = file.path(site_path, fig_url))
    }
    knit_fig_path(source_rmd, md_file, fig_path = fig_path, cache_path = cache_path, quiet = TRUE)
  }
  invisible()
}
