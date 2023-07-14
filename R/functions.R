hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x
}
