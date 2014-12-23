\name{plot.parallelogram}
\alias{plot.parallelogram}
\title{Plot a Parallelogram Object}

\description{ Plot inforce premium levels and rate changes in a
  parallelogram object }

\usage{
\method{plot}{parallelogram}(x, xrange=NULL, fill.colors=c("#7C6C5C", "#583008"),
     text.color="#06401B", written.prem.color="#583008", ...)
}

\arguments{
  \item{x}{A parallelogram object, as produced by the
    \code{\link{Parallelogram}} function}

\item{xrange}{A vector of length two indicating the minimum and maximum
  times to show on the graph.  The default is approximately the earliest
  and latest range change.}

\item{fill.colors}{The alternating colors to use when shading in the
  regions corresponding to each rate level}

\item{text.color}{The color of text to use for each rate level}
\item{written.prem.color}{The color of the line that shows the
  written premium level}

\item{...}{Placeholder required by R}
}

\value{
  Returns a \code{\link{ggplot2}} plot object showing the rate level of
  inforce premium.  The rate of premium written will be shown as a
  dotted line.
}

\seealso{
  \code{\link{Parallelogram}},
  \code{\link{WrittenByPeriod}},
  \code{\link{EarnedByPeriod}},
  \code{\link{ggplot}},
}

\examples{


##### Example 1: constant rate of written premium (classic parallelogram shape)

rate.change.df <- data.frame(year=c(2002, 2003.5, 2004.1, 2004.4, 2005.5),
                             rate.change=c(.07, -.03, .12, .02, .1))
term.len <- 1
p <- Parallelogram(rate.change.df, term.len)
plot(p)

periods.out <- 2001:2006
WrittenByPeriod(p, periods.out)
EarnedByPeriod(p, periods.out)


##### Example 2: Same rate changes as Ex1, but variable rate of written premium

# variable written premium
written.df <- data.frame(year=2001:2005,
                         written=c(30, 45, 75, 30, 55))
p <- Parallelogram(rate.change.df, term.len, written.df)
plot(p)

EarnedByPeriod(p, periods.out)
WrittenByPeriod(p, periods.out)

}
