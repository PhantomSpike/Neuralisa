function db = get_db(data)

xrms = rms(data); %Calculate rms of the signal
Q = 2*(10^-5); %The constant to which all sound levels are compared ~ 20uPa to get the result in sound pressure level (SPL)
db = 20*log10(xrms/Q); %Get the dB