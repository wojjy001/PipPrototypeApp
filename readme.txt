Packages required for this application:
	- shiny
	- shinydashboard
	- ggplot2
	- grid
	- plyr
	- compiler
	- rmarkdown

Need to alter 2 directories in the server.R file

dir = directory where the application is saved

pandocdir = directory where pandoc is saved (a
package required for writing to a word document -
which is automatically installed when "rmarkdown" is 
installed.  If using Windows and have installed RStudio
as normal, then you shouldn't need to change pandocdir.