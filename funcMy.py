def find_nearest(array, value):
    import numpy as np
    idx = (np.abs(array-value)).argmin()
    return idx

def damper(tSoft, time, dtype="exp"):
    import math
    x = time / tSoft
    if dtype == "exp":
        coeff = (math.exp(x**2)-1)/(math.exp(1)-1)
    elif dtype == "exp2":
        coeff = 1.0 / (1 + math.exp(-6*x))
    elif dtype == "linear":
        coeff = x
    else:
        coeff = 1
    return coeff

def rotate2D(origin, point, angle):
    import numpy as np
    angle = np.radians(angle[0])
    ox, oy = origin
    px, py = point
    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy

def rotate3D(origin, point, angle):
    import numpy as np
    alpha = np.radians(angle[0])
    beta = np.radians(angle[1])
    gamma = np.radians(angle[2])
    ox, oy, oz = origin
    px, py, pz = point
    nx = ox + \
         (px - ox) * (np.cos(beta) * np.cos(gamma)) + \
         (py - oy) * (np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma)) + \
         (pz - oz) * (np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma))

    ny = oy + \
         (px - ox) * (np.cos(beta) * np.sin(gamma)) + \
         (py - oy) * (np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma)) + \
         (pz - oz) * (np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma))

    nz = oz + \
         (px - ox) * (-np.sin(beta)) + \
         (py - oy) * (np.sin(alpha) * np.cos(gamma)) + \
         (pz - oz) * (np.cos(alpha) * np.cos(beta))
    return nx, ny, nz

def rotate3Dv2(origin, point, angle):
    import numpy as np
    alpha = np.radians(angle[0])
    beta = np.radians(angle[1])
    gamma = np.radians(angle[2])
    ox, oy, oz = origin
    px, py, pz = point
    nx = ox + \
         (px - ox) * (np.cos(alpha) * np.cos(beta)) + \
         (py - oy) * (np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)) + \
         (pz - oz) * (np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma))

    ny = oy + \
         (px - ox) * (np.sin(alpha) * np.cos(beta)) + \
         (py - oy) * (np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma)) + \
         (pz - oz) * (np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma))

    nz = oz + \
         (px - ox) * (-np.sin(beta)) + \
         (py - oy) * (np.cos(beta) * np.sin(gamma)) + \
         (pz - oz) * (np.cos(beta) * np.cos(gamma))
    return nx, ny, nz

def rao(beta_path, scale, wave_height, freq, angles_type="degrees"):
    import math
    import numpy as np
    gravity = 9.81
    data = np.loadtxt(beta_path).transpose()
    waveA = wave_height / 2
    waveL = (wave_height / 0.17) ** (4 / 3)
    waveS = wave_height / waveL
    ind = np.where(data[0] == freq)[0][0]

    omegaN = data[0][ind]
    periodN = 2 * math.pi / omegaN
    periodM = periodN * math.sqrt(scale)
    omegaM = 2 * math.pi / periodM
    wnN = math.pow(omegaN, 2) / gravity
    wnM = math.pow(omegaM, 2) / gravity
    cN = omegaN / wnN
    cM = omegaM / wnM

    # RAO initilisation
    xi = data[1][ind] * waveA
    eta = data[2][ind] * waveA
    zeta = data[3][ind] * waveA
    theta = data[4][ind] * omegaM ** 2 / gravity * waveA
    psi = data[5][ind] * omegaM ** 2 / gravity * waveA
    chi = data[6][ind] * omegaM ** 2 / gravity * waveA
    xiPhase = data[7][ind]
    etaPhase = data[8][ind]
    zetaPhase = data[9][ind]
    thetaPhase = data[10][ind]
    psiPhase = data[11][ind]
    chiPhase = data[12][ind]

    if angles_type == "degrees":
        theta = np.degrees(theta)
        psi = np.degrees(psi)
        chi = np.degrees(chi)

    print("Waves motion")
    print(f"Phase Velocity:\n\tNatural = {cN:.4f} (m/s)\n\tModel = {cM:.4f} (m/s)")
    print(f"Frequency:\n\tNatural = {omegaN:.4f} (rad/s)\n\tModel = {omegaM:.4f} (rad/s)")
    print(f"Period:\n\tNatural = {periodN:.4f} (s)\n\tModel = {periodM:.4f} (s)")

    print("Linear motions (meters)")
    print("xi \t|\t eta    \t|\t zeta")
    print(f"{xi:.4f} \t|\t {eta:.4f} \t|\t {zeta:.4f}\n")

    print(f"Angular motions ({angles_type})")
    print("theta \t|\t psi    \t|\t chi")
    print(f"{theta:.4f} \t|\t {psi:.4f} \t|\t {chi:.4f}\n")

    return cM, periodM, omegaM, \
           xi, eta, zeta, \
           theta, psi, chi, \
           xiPhase, etaPhase, zetaPhase, \
           thetaPhase, psiPhase, chiPhase

def experiment(beta_path, scale, linear_amplitude, angular_amplitude, freq):
    import math
    import numpy as np
    gravity = 9.81
    data = np.loadtxt(beta_path).transpose()
    ind = np.where(data[0] == freq)[0][0]

    omegaN = data[0][ind]
    periodN = 2 * math.pi / omegaN
    periodM = periodN * math.sqrt(scale)
    omegaM = 2 * math.pi / periodM
    wnN = math.pow(omegaN, 2) / gravity
    wnM = math.pow(omegaM, 2) / gravity
    cN = omegaN / wnN
    cM = omegaM / wnM

    # RAO initilisation
    xi = data[1][ind] * linear_amplitude
    eta = data[2][ind] * linear_amplitude
    zeta = data[3][ind] * linear_amplitude
    theta = data[4][ind] * angular_amplitude
    psi = data[5][ind] * angular_amplitude
    chi = data[6][ind] * angular_amplitude
    xiPhase = data[7][ind]
    etaPhase = data[8][ind]
    zetaPhase = data[9][ind]
    thetaPhase = data[10][ind]
    psiPhase = data[11][ind]
    chiPhase = data[12][ind]

    print("Waves motion")
    print(f"Phase Velocity:\n\tNatural = {cN:.4f} (m/s)\n\tModel = {cM:.4f} (m/s)")
    print(f"Frequency:\n\tNatural = {omegaN:.4f} (rad/s)\n\tModel = {omegaM:.4f} (rad/s)")
    print(f"Period:\n\tNatural = {periodN:.4f} (s)\n\tModel = {periodM:.4f} (s)")

    print("Linear motions (meters)")
    print("xi \t|\t eta    \t|\t zeta")
    print(f"{xi:.4f} \t|\t {eta:.4f} \t|\t {zeta:.4f}\n")

    print("Angular motions (degrees)")
    print("theta \t|\t psi    \t|\t chi")
    print(f"{theta:.4f} \t|\t {psi:.4f} \t|\t {chi:.4f}\n")

    return cM, periodM, omegaM, \
           xi, eta, zeta, \
           theta, psi, chi, \
           xiPhase, etaPhase, zetaPhase, \
           thetaPhase, psiPhase, chiPhase

def rotation6DoF(time, transAmp, transOmega, transPhase, rotAmp, rotOmega, rotPhase):
    import numpy as np
    transX = transAmp[0] * np.sin(transOmega[0] * time + transPhase[0])
    transY = transAmp[1] * np.sin(transOmega[1] * time + transPhase[1])
    transZ = transAmp[2] * np.sin(transOmega[2] * time + transPhase[2])

    rotX = rotAmp[0] * np.sin(rotOmega[0] * time + rotPhase[0])
    rotY = rotAmp[1] * np.sin(rotOmega[1] * time + rotPhase[1])
    rotZ = rotAmp[2] * np.sin(rotOmega[2] * time + rotPhase[2])
    return transX, transY, transZ, rotX, rotY, rotZ


def newSignal(x, y):
    from scipy import signal, interpolate
    from numpy import diff, min, arange
    y = signal.detrend(y)
    dxnew = min(diff(x))
    # dxnew = (x[1] - x[0])
    xnew = arange(x[0], x[len(x) - 1], dxnew)
    f = interpolate.interp1d(x, y, kind='linear')
    ynew = f(xnew)
    return xnew, ynew

def fftAnalyze(x, y):
    from scipy import fft
    from numpy import abs, argmax, ones
    y_fft = fft.fft(y)
    dx = x[1] - x[0]
    freqs = fft.fftfreq(len(x), d=dx)
    #freqs = np.linspace(0, 1.0 / (2.0 * dx), len(x) // 2)
    x_fft = freqs[:len(x) // 2]
    y_fft = 2.0 / len(x) * abs(y_fft[:len(x) // 2])
    y_fft2 = abs(fft.fft(y)/len(x))
    sigFreq = freqs[argmax(y_fft)]
    sigAmp = 2*y_fft[argmax(y_fft)]
    print("Significant Amplitude = {}, Frequency = {}".format(sigAmp, sigFreq))
    meanAmpl = ones(len(x)) * sigAmp
    return x_fft[1:], y_fft[1:], meanAmpl[1:], sigFreq, sigAmp
    #return x_fft, y_fft, meanAmpl, sigFreq, sigAmp

def fftAnalyze2(x, y, nperseg=512, overlap=128):
    from scipy import signal
    from numpy import argmax
    from math import sqrt
    f, Pxx_den = signal.welch(y, fs=1/(x[1]-x[0]), nperseg=1024, scaling='spectrum')
    sigFreq = f[argmax(Pxx_den)]
    sigAmp = 2*sqrt(Pxx_den[argmax(Pxx_den)])
    print("Significant Amplitude = {}, Frequency = {}".format(sigAmp, sigFreq))