---
layout: post
title: Markdown Cheatsheet  
published: true
excerpt: >
  A cheatsheet and showcase of Markdown syntax meant to show how different 
  elements will be rendered.  
category: web
tags: r markdown
---



This website is generated from [Markdown](https://en.wikipedia.org/wiki/Markdown) files using the static site generator [Jekyll](jekyllrb.com). When editing the templates and CSS that define how the Markdown files are converted to html, I use this document to see how all the various elements will appear in the browser. Pieces of this were adapted from [this cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

# Paragraphs

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Cras accumsan, libero nec facilisis rhoncus, velit turpis vestibulum orci, a pharetra nisi libero vel orci. Vivamus ante lacus, sagittis id est nec, rutrum venenatis dui. In vel elit dapibus leo posuere volutpat. Integer eleifend risus ut ligula ultricies, sit amet accumsan quam lobortis. Aliquam et rhoncus quam, eu scelerisque nulla.

Nunc magna nisi, rutrum sit amet mollis in, maximus ac erat. Suspendisse fringilla sit amet metus at suscipit. Maecenas ultricies venenatis metus et feugiat. Sed lobortis congue tristique. Nulla pharetra egestas feugiat. Curabitur rhoncus eu arcu vel volutpat. Integer cursus lorem nunc, nec fermentum ex consectetur auctor. Nullam a tortor nec diam suscipit porta sed sit amet felis.

Donec congue, nisi id tincidunt finibus, sem risus accumsan quam, et cursus diam tellus eget lorem. Integer aliquam bibendum mi in tincidunt. Vestibulum sit amet convallis nibh. Phasellus rutrum est sed urna lacinia varius. Pellentesque lobortis porttitor turpis, sit amet porttitor mauris aliquam ac. Donec scelerisque hendrerit velit dignissim tincidunt. Donec varius pharetra magna, sit amet consequat ex condimentum vel. Nulla sapien mauris, fringilla eget lectus sit amet, porttitor ultrices est. Nullam justo metus, viverra at tortor tempus, malesuada molestie nulla. Aliquam mattis a nisl et elementum.

# Headings

# Heading 1

## Heading 2

### Heading 3

#### Heading 4

##### Heading 5

###### Heading 6

# Styling Text

Emphasis, aka italics, with *asterisks* or _underscores_.

Strong emphasis, aka bold, with **asterisks** or __underscores__.

Combined emphasis with **asterisks and _underscores_**.

# Lists

## Ordered

My favourite foods:

1. Mapo tofu
2. Spanakopite
3. Poutine
4. Donec congue, nisi id tincidunt finibus, sem risus accumsan quam, et cursus diam tellus eget lorem. Integer aliquam bibendum mi in tincidunt. Vestibulum sit amet convallis nibh.
5. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Cras accumsan, libero 
   nec facilisis rhoncus, velit turpis vestibulum orci, a pharetra nisi libero vel orci.

## Unordered

A list of names:

- Matt
- Josh
- John
- Sarah
- Katie
- Phasellus rutrum est sed urna lacinia varius. Pellentesque lobortis porttitor turpis, sit amet porttitor mauris aliquam ac. Donec scelerisque hendrerit velit dignissim tincidunt. Donec varius pharetra magna, sit amet consequat ex condimentum vel. Nulla sapien mauris, fringilla eget lectus sit amet, porttitor ultrices est. Nullam justo metus, viverra at tortor tempus, malesuada molestie nulla. Aliquam mattis a nisl et elementum.

## Nested lists

1. First item
    1. First sub-item
    2. Second sub-item
2. Second item
    - Something
    - Soemthing else
3. Third item

- Some item
    - Sub point
    - Another sub-point
- Another item
    1. First
    2. Second
- Last thing

# Links

This is a link to the [wikipedia page on Markdown](https://en.wikipedia.org/wiki/Markdown). And here's a page about [GitHub flavoured Markdown](https://help.github.com/articles/github-flavored-markdown/). Or maybe you just want to visit [Google](https://www.google.com).

# Code

Inline snippets of `code` are placed between backticks. Code chunks can be fenced in triple backticks, which will use a fixed width font, preserve spaces, and add syntax highlighting:

```html
<html>
  <head>
    <title>Markdown Cheatsheet</title>
  </head>
</html>
```

```r
# This is a really really really long comment to see what happens when lines of code overflow box, will it wrap of scroll???
s <- "Hello World!"
print(s)
```

## Dynamic Code Chunks

`knitr` will run code chunks dynamically in R and include both the source code and the output in the markdown document.


```r
x <- runif(100)
y <- 2 * x + 1 + rnorm(100, mean = 0.5, sd = 0.1)
plot(x, y)
```

<img src="/figures//markdown-cheatsheet_unnamed-chunk-1-1.svg" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />

Inline code can also be dynamic, for example this is the mean miles per gallon of cars in the `mtcars` dataset: 20.091.

# Mathematical Expressions

Mathematical expressions can be rendered in the browswer with [MathJax](https://www.mathjax.org/). To define a display equation in Markdown it should be surrounded by double dollar signs (`$$`):

$$
\oint_{\partial \Sigma} \mathbf{B} \cdot \mathrm{d}\boldsymbol{\ell} = \mu_0 \iint_{\Sigma} \mathbf{J} \cdot \mathrm{d}\mathbf{S} + \mu_0 \varepsilon_0 \frac{\mathrm{d}}{\mathrm{d}t} \iint_{\Sigma} \mathbf{E} \cdot \mathrm{d}\mathbf{S} 
$$

Inline equations are surrounded with single dollar signs, or using the Redcarpet or kramdown markdown renderer in Jekyll, `\\( E = m c^2 \\)` yields \\( E = m c^2 \\).

# Images

A 400x400 image:

![400x400 Placeholder](http://placehold.it/400)

A smaller, rectangular (200x50) image:

![smaller image](http://placehold.it/350x150)

Finally, a huge (1000x500) image:

![huge image](http://placehold.it/1000x500)

# Table

A simple table:

First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  | Content Cell

More complicated table with numbers and alignment.

| Left-Aligned  | Center Aligned  | Right Aligned |
| :------------ |:---------------:| -----:|
| col 3 is      | some wordy text | $1600 |
| col 2 is      | centered        |   $12 |
| zebra stripes | are neat        |    $1 |

Finally, a table generated by R with `knitr::kable()`:


```r
knitr::kable(mtcars[1:10, 1:6], digits = 1)
```



|                  |  mpg| cyl|  disp|  hp| drat|  wt|
|:-----------------|----:|---:|-----:|---:|----:|---:|
|Mazda RX4         | 21.0|   6| 160.0| 110|  3.9| 2.6|
|Mazda RX4 Wag     | 21.0|   6| 160.0| 110|  3.9| 2.9|
|Datsun 710        | 22.8|   4| 108.0|  93|  3.8| 2.3|
|Hornet 4 Drive    | 21.4|   6| 258.0| 110|  3.1| 3.2|
|Hornet Sportabout | 18.7|   8| 360.0| 175|  3.1| 3.4|
|Valiant           | 18.1|   6| 225.0| 105|  2.8| 3.5|
|Duster 360        | 14.3|   8| 360.0| 245|  3.2| 3.6|
|Merc 240D         | 24.4|   4| 146.7|  62|  3.7| 3.2|
|Merc 230          | 22.8|   4| 140.8|  95|  3.9| 3.1|
|Merc 280          | 19.2|   6| 167.6| 123|  3.9| 3.4|

# Horizontal Rule

A horizontal rule can be inserted with three hyphens,

---

asterisks,

***

or underscores

___

# Blockquotes

A couple classic evolution quotes. One from Darwin:

> Thus, from the war of nature, from famine and death, the most exalted object which we are capable of conceiving, namely, the production of the higher animals, directly follows. There is grandeur in this view of life, with its several powers, having been originally breathed into a few forms or into one; and that, whilst this planet has gone cycling on according to the fixed law of gravity, from so simple a beginning endless forms most beautiful and most wonderful have been, and are being, evolved.

And one from Wallace:

> Every species has come into existence coincident both in space and time with a
> pre-existing closely allied species.
