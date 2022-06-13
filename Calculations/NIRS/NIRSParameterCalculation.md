# NIRS Parameter Calculation Algorithm & Explanation.

Below is a summary of the implementation of the algorithm in both Matlab and JavaScript. Each step has both Matlab and JavaScript codes to be able to compare them.

The JavaScript code is an exact implementation of the Matlab code but as a class (object oriented programming).

---

---

### Step 1 - Coefficient Calculation and Variable Definition.

Here we define the initial coefficients and variables used in the heomodynamics and TOI calculation.

##### **MATLAB Code**

Load coefficients for hemoglobin and TOI calculation. Matrix`c_beta` is only used for TOI calculation.

```
load e_coef
c_beta = [-1.2132   -0.0728    1.8103    1.1433  -11.5816];
```

##### **JS Code**

As JavaScript does have native support for matrices, the code is implemented using multi-dimensional arrays. Here `this` refers to the instance of the class that does the NIRS parameter calculation. Variables are defined globally within the class to be able to use them inside all methods of the class.

```
this.wavelengths = new Uint16Array([950, 730, 810, 850, 650]);
this.c_beta = new Float32Array([-1.2132, -0.0728, 1.8103, 1.1433, -11.5816]);
this.ADCRes = 12;
this.ADCMaxVal = Math.pow(2, this.ADCRes) - 1;
this.DACRes = 8;
this.DACMaxVal = Math.pow(2, this.DACRes) - 1;
```

##### **MATLAB Code**

Define LED wavelengths and calculate the wave index used to grab the appropriate coefficients from the `e_coef` matrix of coefficients.

```
wavelength = [950 730 810 850 650]; // V5's wavelengths in the right order
waveidx = round((wavelength - 599.9)/.1);
```

Calculate the coefficient matrix

```
HBcoef(:,1) = e_coef(waveidx,3);
HBcoef(:,2) = e_coef(waveidx,2);
A = -pinv(HBcoef);
```

##### **JS Code**

Calculate the wave index.

```
this.waveIndex = this.wavelengths.map((wavelength) => Math.round((wavelength - 599.9) / 0.1));
```

Calculate the HB Coef. Here `col1`, `col2`, `col3` refers to the columns of the `e_coef` matrix defined in the MATLAB implementation. In JavaScript, `e_coef` has been split into 3 Float32 arrays named `col1`, `col2`, and `col3`.

```
const HBCoef1 = new Float32Array(this.waveIndex.length);
const HBCoef2 = new Float32Array(this.waveIndex.length);

// Loop through each wave index and find the coefficients from the e_coef matrix.
this.waveIndex.forEach((waveIdx, i) => {
    HBCoef1[i] = col3[waveIdx];
    HBCoef2[i] = col2[waveIdx];
});

this.HBCoef = [HBCoef1, HBCoef2];
```

Calculate A.

```
const _matrixA: number[][] = []; // Matrix A is a 2d matrix of numbers

// Creates a 5 x 2 matrix to be used to get its pinv
this.HBCoef[0].forEach((_coef, i) => {
    _matrixA.push([this.HBCoef[0][i], this.HBCoef[1][i]]);
 });

// Since JS does not pseudo matrix inverse, we use a widely used library called ML matrix
// to bring matrix calculation to JS. Here we define a Matrix in the ML matrix library based
// on their API.
const matrixA = new Matrix(_matrixA);

// The library calculates the pseudo inverse of the matrix, does a matrix multiplication by -1
// And will return the data as a JavaScript 2d array using the `.data` getter.
this.A = pseudoInverse(matrixA).multiply(-1).data;
```

---

### Step 2 - Calculating Hemoglobin Relative Concentration Values.

##### **MATLAB Code**

```
% Calculating hemoglobin concentration values
% W1, ..., W5 are raw wavenlegth values reported by the sensor at each time
% instant (each data point). B is the baseline value. If baseline removal is not needed, set
% B=0.
B = 0;
W = nirsData(1:5, :)';
data = (W-B)/4096; % 4096 is the 12 bit ADC max value.
data = max(data, 0.001)';
HBconc = A*log(data);
O2Hb = HBconc(1,:)';
HHb  = HBconc(2,:)';
tHb  = O2Hb + HHb;
```

##### **JS Code**

Here sample is a number array of 6 elements with the following format: `[ambient, w1, w2, w3, w4, w5]`

```
// Subtract the baseline from each raw value and divide by (4096 - ADC)
// The first element of the sample array is the ambient/baseline
for (let i = 1; i < this.PDChannels; i += 1) {
    sample[i] = (sample[i] - sample[0]) / this.ADCMaxVal;

    // If the value is less than 0.001, replace it with 0.001
    if (sample[i] < 0.001) sample[i] = 0.001;

    sample[i] = Math.log(sample[i]);
}

let O2Hb = 0;
let HHb = 0;

// Use matrix dot product to calculate O2HB and HHb
for (let i = 1; i < this.PDChannels; i++) {
    O2Hb += sample[i] * this.A[0][i - 1];
    HHb += sample[i] * this.A[1][i - 1];
}

// Final values
O2Hb = Math.abs(O2Hb);
HHb = Math.abs(HHb);
```

---

### Step 3 - Calculate the Matrix of Coefficient for TOI.

##### **MATLAB Code**

The code below calculate the matrix of coefficients for TOI calculation. Here `led` is the driving intensity of each LED, `255` is the 8 DAC max value, and `4096` is the 12 bit ADC max value.

```
led = double(led);
Amps_coef = [led(1), led(2), led(4), led(5)]/255*4096; % Replace 255*4096 with 127*16383 for Teliatry sensor
L_coefreq = [W(:,1), W(:,2), W(:,4), W(:,5)] ./ W(:,3);
```

##### **JS Code**

Here sample is a number array of 6 elements with the following format: `[ambient, w1, w2, w3, w4, w5]`

```
// Remove 3rd element of the array - Used to normalize
const Amp_coef = new Float32Array([
    (this.LEDIntensities[0] / this.DACMaxVal) * this.ADCMaxVal,
    (this.LEDIntensities[1] / this.DACMaxVal) * this.ADCMaxVal,
    (this.LEDIntensities[3] / this.DACMaxVal) * this.ADCMaxVal,
    (this.LEDIntensities[4] / this.DACMaxVal) * this.ADCMaxVal,
]);

// Normalize based on one Raw PD (ADC) value - LED 3 Wavelength chosen here
const L_coefreq = new Float32Array([
    sample[1] / sample[3],
    sample[2] / sample[3],
    sample[4] / sample[3],
    sample[5] / sample[3],
]);
```

---

### Step 4 - Calculate TOI.

Uses the coefficients calculated in step 3 and the `c_beta` matrix to calculate the TOI.

```
OD_test = -log(L_coefreq./Amps_coef)';
TOI = c_beta*[OD_test; ones(1, size(OD_test,2))];
TOI = abs(TOI')*100;
```

##### **JS Code**

```
// Create an array of size: number of LEDs and fill it with 1.0
// This is used to make the array the same size as the number of LEDs
// while the L_coefreq and Amp_coef size is 1 less than the number of LEDs
// fill the rest with 1.0
const OD_Test = new Float32Array(this.NUM_OF_LEDs).fill(1.0);
let TOI = 0;

for (let i = 0; i < this.NUM_OF_LEDs; i++) {
    // Only do the operation for the same size(4) of L_coefreq and Amp_coef arrays
    if (i !== this.NUM_OF_LEDs - 1) OD_Test[i] = Math.log(L_coefreq[i] / Amp_coef[i]) * -1;

    TOI += this.c_beta[i] * OD_Test[i];
}

// Take the absolute value and multiply by 100 to get a positive percentage
TOI = TOI * 100; // Absolute value removed
```

---
