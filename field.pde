// GLOBAL VARIABLES
// ===========================================================================
/* Diff. of Space */ float amplitude;
/* Cell Number in a Row */ int number_of_cells;
/* Diff. of Time */ float dt;
/* Pix. Number W/H */ int pixels;
/* Table of Results */ Table results;
/* Maximum number of time steps, if used. */ int time_steps;
/* Maximum density for triggers */ float maximum_rho;
/* Decl. Grid */ cell[][] field;
/* Trigger Array */ ArrayList<trigger> density;
/* Const. Value for Init. Conditions */ float init_conditions;
/* Very 1st (Test) Trigger */ trigger new_trigger;
/* Time Randomness */ int Frequency;
/* Vector to classify triggers according to sizes */ int[] classify_triggers;
/* Bool Variable for identifying trigger and cells */ boolean trigger_or_cell;
/* Determined trigger dimensions */ int trigger_size;
/* Time step counter */ int time;
/* Max. temperature in heat field */ float temp_max;
/* Number of subdivisions for the data.csv file */ int precision;
// ===========================================================================
void setup() {
	// SETTING UP THE SYSTEM . . . 
	// Time aspects
		time = 0;
		Frequency = 2;
		dt = 0.1;
		time_steps = 1;
	// Space aspects
		number_of_cells = 200;
		pixels = 200;
		size(pixels,pixels);
		amplitude = pixels / number_of_cells;
		field = new cell[number_of_cells][number_of_cells];
	// Density aspects
		maximum_rho = 1;
		trigger_size = 2 * number_of_cells / 100;
		init_conditions = 254;
		trigger_or_cell = true; // Trigger -> true, Cell -> false
		density = new ArrayList<trigger>();
		temp_max = ((float)Math.pow(((2 * trigger_size) + 1),2) * init_conditions) + 1;
	// Data aspects
		precision = 25;
		classify_triggers = new int[precision];
		for (int i = 0; i < precision; i++) classify_triggers[i] = 0;
		// Setting Table for Collecting Data.
		results = new Table();
		results.addColumn("Time"); // Time step number (t)
		for (int j = 0; j < precision; j++) results.addColumn("Size: " + temp_max / precision * (j + 1));
		 // Evaluating trigger size (p) [interval]
	// CONFIGURING FIELD . . . 
	for (int i = 0;  i < number_of_cells; i++) for (int j = 0; j < number_of_cells; j++) {
		field[i][j] = new cell(i * amplitude, j * amplitude, init_conditions, i, j);
	}
	/* 1st Trigger Formed */ form_trigger(); 
}
// FUNCTION FOR CREATING TRIGGER EVENTS 
void form_trigger() {
	int i = (int)random(number_of_cells - 2) + 2;
	int j = (int)random(number_of_cells - 2) + 2;
	trigger new_trigger = new trigger(i * amplitude,j * amplitude,random(maximum_rho),i,j,trigger_size);
	density.add(new_trigger);
	new_trigger.trigger_region();
}
// FUNCTION FOR HEAT DIFFUSION SIMULATION
float heat_equation(int i, int j) {
	float new_temperature = (dt / (amplitude * amplitude)) * 
							(field[i + 1][j].temperature + 
		 					 field[i - 1][j].temperature + 
		 					 field[i][j + 1].temperature + 
		 					 field[i][j - 1].temperature - 
		 					 (4 * field[i][j].temperature));
	return new_temperature;
}
void draw() {
	if (time < pow(10,10)) { /* Setting time limitation */
		// Boundary conditions (constant value for rows and columns i,j = {0,n - 1})
		for (int i = 0;  i < number_of_cells; i++) {
			field[i][0].temperature = 0;
			field[i][number_of_cells - 1].temperature = 0;
			field[number_of_cells - 1][i].temperature = 0;
			field[0][i].temperature = 0;
		}
		if ((int)random(Frequency) == 1) form_trigger(); // Placing triggers randomly in time ...
		// Applying Heat Equation
			for (int i = 1; i < number_of_cells - 1; i++) for (int j = 1; j < number_of_cells - 1; j++) {
				for (int k = 0; k < density.size(); k++) {
					if ((density.get(k).indexi == i && density.get(k).indexj == j)      || 
					    (density.get(k).indexi == i + 1 && density.get(k).indexj == j)  ||
					    (density.get(k).indexi == i - 1 && density.get(k).indexj == j)  || // TO ENSURE NOT COMPUTING HEAT DIFFUSION
					    (density.get(k).indexi == i && density.get(k).indexj == j + 1)  || // WHEN THERE IS A TRIGGER
					    (density.get(k).indexi == i && density.get(k).indexj == j - 1)) {
						trigger_or_cell = true;
						break; // If there is a trigger, we do not simulate heat diffusion.
					}
					else trigger_or_cell = false;
				}
				if (trigger_or_cell == false) field[i][j].temperature += heat_equation(i,j); // SIMULATE HEAT DIFFUSION
			}
		// Printing method for all the field 
		for (int i = 0;  i < number_of_cells; i++) for (int j = 0; j < number_of_cells; j++) field[i][j].exist(); 
		} 
		time++;
		// DATA ORGANIZATION
		/* Every 100 time steps */ if (time % 100 == 0) { 
			for (int i = 0; i < precision; i++) classify_triggers[i] = 0; // Setting values to 0 to restart summing trigger number.
			for (int i = 0; i < density.size(); i++) for (int j = 0; j < precision; j++) {
				// WE CONSIDER A TRIGGER BEING WITHIN AN INTERVAL WHICH IS DETERMINED BY:
					// * precision
					// * RHO of the trigger
				if (density.get(i).rho >= j * temp_max / precision && density.get(i).rho <= (j + 1) * temp_max / precision) {
					classify_triggers[j]++;
				}
			}
			TableRow newRow = results.addRow();
			newRow.setInt("Time",time);
			// A ROW FOR EACH INTERVAL IS GENERATED
			for (int j = 0; j < precision; j++) {
				newRow.setInt("Size: " + temp_max / precision * (j + 1),classify_triggers[j]);
			}
			// DATA SAVED IN "data.csv" file
			saveTable(results,"data.csv");
			print(time + "," + density.size() + "\n"); // TERMINAL NOTIFICATION 
		}
}
// ~ C E L L  C L A S S ~ 
// ====================== 
class cell {
	float posx, posy, temperature;
	int indexi, indexj;
	cell(float posxa, float posya, float temperaturea, int indexia, int indexja) {
		posx = posxa;
		posy = posya;
		temperature = temperaturea;
		indexi = indexia;
		indexj = indexja;
	}
	void exist() {
		float temp_max = ((float)Math.pow(((2 * trigger_size) + 1),2) * init_conditions) + 1;
		fill(0,0,temperature * 0.95 + 20);
		stroke(0,0,temperature * 0.95 + 20);
		rect(posx,posy,posx + amplitude, posy + amplitude);
	}
};
// ~ T R I G G E R  C L A S S ~
// ============================
class trigger {
	float posx, posy, rho;
	int dimensions;
	int indexi, indexj;
	trigger(float posxa, float posya, float rhoa, int indexia, int indexja, int dimensionsa) {
		posx = posxa;
		posy = posya;
		rho = rhoa;
		indexi = indexia;
		indexj = indexja;
		dimensions = dimensionsa;
	}
	void trigger_region() {
		float const_rho = rho;
		float Probability = 0;
		float ideal_value = 0;
		// Covering "Trigger Area" to direct temperature to 
		for (int i = indexi - dimensions; i <= indexi + dimensions; i++) for (int j = indexj - dimensions; j <= indexj + dimensions; j++) {
			if (i > 1 && i < number_of_cells - 1 && j > 1 && j < number_of_cells - 1) if (i != indexi || j != indexj) {
				// % OF RHO DETERMINES HOW MUCH TEMPERATURE IS SUCKED BY THE TRIGGER
				rho += const_rho / maximum_rho * field[i][j].temperature;
				ideal_value += const_rho / maximum_rho * init_conditions;
				// FIELD TEMPERATURE IS REDUCED
				field[i][j].temperature = (maximum_rho - const_rho) / maximum_rho * field[i][j].temperature;
			}
		}
		Probability = rho / ideal_value;
		if (Probability > 0.5) field[indexi][indexj].temperature = rho / (((float)Math.pow(((2 * dimensions) + 1),2) * init_conditions) + 1) * 255;
		else field[indexi][indexj].temperature = -pow(10,10);
	}
};
