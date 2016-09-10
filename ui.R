#Define UI for PipPrototypeApp2
library(shinydashboard)
#-------------------------------------------------------------------------------------
#Header
header <- dashboardHeader(
	title = "Piperacillin Application",
	titleWidth = 250
)  #Brackets closing "header"

#Sidebar
sidebar <- dashboardSidebar(
	width = 250,  #Width of sidebar
	sidebarMenu(
		menuItem("1. Patient Information", tabName = "patient", icon = icon("child")),
		menuItem("2. Dosing Information", tabName = "dosing", icon = icon("medkit")),
		menuItem("3. Plot and Numerical Output", tabName = "graphs", icon = icon("line-chart")),
		menuItem("About", tabName = "about", icon = icon("question"),
			menuSubItem("Model", tabName = "model", icon = icon("angle-right"))
		)
	)	#Brackets closing "sidebarMenu"
)	#Brackets closing "dashboardSidebar"

#Body
body <- dashboardBody(
	tabItems(
		tabItem(tabName = "patient",
			fixedRow(
				column(4,
					textInput("PHARMI", "Pharmacist's Initials:", "Your Initials")
				)
			),	#Brackets closing "fixedRow"
			h4(strong("Patient Details:")),
			fixedRow(
				column(4,
					textInput("FNAME", "First Name:", "First Name"),
					textInput("LNAME", "Last Name:", "Last Name")
				),
				column(4,
					numericInput("URN", "Unit Record Number (URN):", value = 999999, step = 1)
				)
			),	#Brackets closing "fixedRow"
			br(),
			h4(strong("Covariate Information:")),
			fixedRow(
				column(4,
					dateInput("BDATE", "Date of Birth (DD-MM-YYYY):", value = "1930-08-30", format = "dd-mm-yyyy", startview = "year"),
					numericInput("WT", "Weight (kg):", value = 61.1, step = 0.1)
				),
				column(4,
					selectInput("SEX", "Gender:", choices = c("Male" = 1, "Female" = 2), selected = 1),
					numericInput("SECR", "Serum Creatinine (µmol/L):", value = 60.2, step = 0.1),
					textOutput("CrCLText")
				)
			),	#Brackets closing "fixedRow"
			h4(strong("Previous Information:")),
			fixedRow(
				column(8,
					div(style = "overflow-x: auto", tableOutput("prevTable"))
				)
			)  #Brackets closing "fixedRow"
		),
		tabItem(tabName = "dosing",
			h4(strong("Prescribed Dosing Information:")),
			fixedRow(
				column(4,
					dateInput("DDATE", "Date of Sampled Dose Administration:", value = NULL, format = "dd-mm-yyyy", startview = "month")
				)
			),	#Brackets closing "fixedRow"
			fixedRow(
				column(4,
					numericInput("PDOSE", "Dose Amount (mg):", value = 4000, min = 0, max = 5000, step = 100)
				),
				column(4,
					selectInput("PFREQ", "Dose Frequency (hours):", choices = c("4-hourly" = 1, "6-hourly" = 2,"8-hourly" = 3), selected = 2)
				)
			),  #Brackets closing "fixedRow"
			fixedRow(
				column(4,
					numericInput("PINFD", "Prescribed Infusion Duration (hours):", value = 3, min = 0.5, max = 8, step = 0.5)
				),
				column(4,
					selectInput("NPDOSE", "Number of Previous Doses:", choices = c("No previous doses, this is the first dose" = 1, "1" = 2, "2" = 3, "3 or more previous doses (at steady state)" = 4), selected = 3)
				)
			),  #Brackets closing "fixedRow"
			br(),
			h4(strong("Collected Concentration Information:")),
			fixedRow(
				column(4,
					selectInput("NCONC", "Number of samples collected:", choices = c("1" = 1, "2" = 2, "3" = 3, "4" = 4), selected = 2)
				)
			),  #Brackets closing "fixedRow"
			fixedRow(
				column(4,
					numericInput("PCONC1", "1: Concentration (mg/L)", value = 70, min = 0.01, step = 0.001),
					conditionalPanel(condition = "input.NCONC >= 2",
						numericInput("PCONC2", "2: Concentration (mg/L)", value = 20, min = 0.01, step = 0.001)
					),  #Brackets closing "conditionalPanel" for when NCONC >= 2
					conditionalPanel(condition = "input.NCONC >= 3",
						numericInput("PCONC3", "3: Concentration (mg/L)", value = 25, min = 0.01, step = 0.001)
					),  #Brackets closing "conditionalPanel" for when NCONC >= 3
					conditionalPanel(condition = "input.NCONC == 4",
						numericInput("PCONC4", "4: Concentration (mg/L)", value = 10, min = 0.01, step = 0.001)
					)	#Brackets closing "conditionalPanel" for when NCONC == 4
				),
				column(4,
					numericInput("PTIME1", "1: Time after infusion started (hours):", value = 1, step = 0.5),
					conditionalPanel(condition = "input.NCONC >= 2",
						numericInput("PTIME2", "2: Time after infusion started (hours):", value = 5, step = 0.5)
					),  #Brackets closing "conditionalPanel" for when NCONC >= 2
					conditionalPanel(condition = "input.NCONC >= 3",
						numericInput("PTIME3", "3: Time after infusion started (hours):", value = 2, step = 0.5)
					),  #Brackets closing "conditionalPanel" for when NCONC >= 3
					conditionalPanel(condition = "input.NCONC == 4",
						numericInput("PTIME4", "4: Time after infusion started (hours):", value = 3, step = 0.5)
					)  #Brackets closing "conditionalPanel" for when NCONC == 4
				)
			),  #Brackets closing "fixedRow"
			br(),
			h4(strong("Target Information:")),
			fixedRow(
				column(4,
					div(style = "height: 50vh; overflow-y: auto", selectInput("MIC", "Minimum Inhibitory Concentration (MIC):", choices = c("0.25 mg/L" = 1, "0.5 mg/L" = 2, "1 mg/L" = 3, "2 mg/L" = 4, "4 mg/L" = 5, "8 mg/L" = 6, "16 mg/L" = 7, "32 mg/L" = 8, "64 mg/L" = 9), selected = 4)
				))
			)	#Brackets closing "fixedRow"
		),
		tabItem(tabName = "graphs",
			fixedRow(
				column(2, offset = 9,
					actionButton("ADD", label = "Save Patient Data"),
					align = "right"
				),
				column(1,
					conditionalPanel(condition = "input.ADD",
						h5("Saved!", align = "left")
					)
				)
			),	#Brackets closing "fixedRow"
			box(
				fixedRow(
					column(7,
						h4(strong("Piperacillin Concentration-Time Profile")),
						plotOutput("concPlot1", height = 450)
					),
					column(5,
						h4(strong("Predict concentrations for the next 3 doses")),
						sliderInput("SDOSE1", "New Dose Amount (mg):", value = 4000, min = 0, max = 5000, step = 100, width = 400),
						sliderInput("SINFD1", "New Infusion Duration (hours):", value = 3, min = 0.5, max = 8, step = 0.5, width = 400),
						selectInput("SFREQ1", "New Dose Frequency (hours):", choices = c("4-hourly" = 1, "6-hourly" = 2,"8-hourly" = 3), selected = 1, width = 400),
						hr(),
						h4(strong("Time above MIC, % (95% CI)")),
						tableOutput("micTable1"),
						downloadLink("downloadReport", label = h4(strong("Click here to download patient report (.docx)")))
					),
					align = "center"
				),	#Brackets closing "fixedRow"
				width = 12,
				status = "primary"
			),	#Brackets closing "box"
			box(
				title = "Explore a different dosing regimen (for comparison only)",
				fixedRow(
					column(7,
						h4(strong("Piperacillin Concentration-Time Profile")),
						plotOutput("concPlot2", height = 360)
					),
					column(5,
						h4(strong("Predict concentrations for the next 3 doses")),
						sliderInput("SDOSE2", "New Dose Amount (mg):", value = 4000, min = 0, max = 5000, step = 100, width = 400),
						sliderInput("SINFD2", "New Infusion Duration (hours):", value = 3, min = 0.5, max = 8, step = 0.5, width = 400),
						selectInput("SFREQ2", "New Dose Frequency (hours):", choices = c("4-hourly" = 1, "6-hourly" = 2,"8-hourly" = 3), selected = 1, width = 400),
						hr(),
						h4(strong("Time above MIC, % (95% CI)")),
						tableOutput("micTable2")
					),
					align = "center"
				),	#Brackets closing "fixedRow"
				width = 12,
				status = "primary",
				collapsible = TRUE,
				collapsed = TRUE
			)	#Brackets closing "box"
		),
		tabItem(tabName = "model",
			tableOutput("testTable")
		)
	)  #Brackets closing "tabItems"
)  #Brackets closing "body"
#-------------------------------------------------------------------------------------
#User-interface Object
dashboardPage(header, sidebar, body, skin = "blue")
