def turbulence(velocity, dtype="laminar", L="1", Tu="1", nu="1e-6", nutnu="1"):
    import math
    if dtype == "RAS":
        turbulenceModel = "RAS"

        # Turbulence properties
        Cmu = 0.09
        u1 = Tu / 100 * velocity
        nut = nutnu * nu
        tKE = 3 / 2 * math.pow(u1, 2)
        epsilon = Cmu * math.pow(tKE, 2) / (nutnu * nu)
        omega = epsilon / (tKE * Cmu)
        TuL = Cmu * math.pow(L, 3 / 2) / epsilon

        print("###### Turbulent Model")
        print(f"turbulentKE:\t{tKE}")
        print(f"epsilon:\t{epsilon:.6f}")
        print(f"omega:\t{omega:.6f}")
        print(f"TuL:\t{TuL:.6f}")
    else:
        turbulenceModel = "laminar"
        nut = tKE = epsilon = omega = 0
        print("###### Laminar Model")
    return turbulenceModel, nut, tKE, epsilon, omega

def write6DoF(path, startTime, endTime, tSoft, N, \
              transAmp, transOmega, transPhase, \
              rotAmp, rotOmega, rotPhase, damperType="linear"):
    import numpy as np
    from funcMy import damper
    time = np.linspace(startTime, endTime, N)

    transX = transAmp[0] * np.sin(transOmega[0] * time + transPhase[0])
    transY = transAmp[1] * np.sin(transOmega[1] * time + transPhase[1])
    transZ = transAmp[2] * np.sin(transOmega[2] * time + transPhase[2])

    rotX = rotAmp[0] * np.sin(rotOmega[0] * time + rotPhase[0])
    rotY = rotAmp[1] * np.sin(rotOmega[1] * time + rotPhase[1])
    rotZ = rotAmp[2] * np.sin(rotOmega[2] * time + rotPhase[2])

    with open(path, "w") as f:
        f.write("{}\n(\n".format(N))
        for i in range(N):
            if time[i] < tSoft:
                time_ = time[i]
                damper = damper(tSoft, time[i], dtype=damperType)
                transX_ = transX[i] * damper
                transY_ = transY[i] * damper
                transZ_ = transZ[i] * damper
                rotX_ = rotX[i] * damper
                rotY_ = rotY[i] * damper
                rotZ_ = rotZ[i] * damper
            else:
                time_, transX_, transY_, transZ_, rotX_, rotY_, rotZ_ = \
                    time[i], transX[i], transY[i], transZ[i], rotX[i], rotY[i], rotZ[i]
            f.write(f"({time_:.6f} (({transX_:.6f} {transY_:.6f} {transZ_:.6f}) ({rotX_:.6f} {rotY_:.6f} {rotZ_:.6f})))\n")
        f.write(")")
    f.close()

def setup_openfoam(path):
    with open(f"{path}/initial.dat", "w") as f:
        f.write("np={}\n".format(int(numProc)))
        f.write("hpc={}\n".format(hpc))
        f.write("scale={}\n".format(scale))
        f.write("turbulenceModel={}\n".format(turbulenceModel))
        f.write("tKE={}\n".format(tKE))
        f.write("epsilon={}\n".format(epsilon))
        f.write("nut={}\n".format(nut))
        f.write("CoG=\"({} {} {})\"\n".format(CoG[0], CoG[1], CoG[2]))
        f.write("fill={:.4f}\n".format(fillLevel))
        f.write("omegaM={:.4f}\n".format(omegaM))
        f.write("endTime={:.6f}\n".format(endTime))
        f.write("deltaT={:.6f}\n".format(deltaT))
        f.write("writeInterval={:.6f}\n".format(writeInterval))
        f.write("purgeWrite={}\n".format(purgeWrite))
    f.close()