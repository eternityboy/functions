def writeConfig(path, geometry, solver, fill, CoR, \
                 transAmp, transOmega, transPhase, rotAmp, rotOmega, rotPhase):
    import configparser

    config = configparser.ConfigParser()

    config.add_section("Geometry")
    config.set("Geometry", "scale", f"{geometry[0]:.6f}")
    config.set("Geometry", "length", f"{geometry[1]:.6f}")
    config.set("Geometry", "width", f"{geometry[2]:.6f}")
    config.set("Geometry", "height", f"{geometry[3]:.6f}")
    config.set("Geometry", "h_lc", f"{geometry[4]:.6f}")
    config.set("Geometry", "h_uc", f"{geometry[5]:.6f}")
    config.set("Geometry", "theta_lc", f"{geometry[6]:.6f}")
    config.set("Geometry", "theta_uc", f"{geometry[7]:.6f}")

    config.add_section("Solver")
    config.set("Solver", "startTime", f"{solver[0]:.6f}")
    config.set("Solver", "endTime", f"{solver[1]:.6f}")
    config.set("Solver", "deltaT", f"{solver[2]:.6f}")
    config.set("Solver", "writeInterval", f"{solver[3]:d}")
    config.set("Solver", "cfl", f"{solver[4]:.6f}")
    config.set("Solver", "adaptiveTimeStep", f"{solver[5]}")
    config.set("Solver", "writeFrequency", f"{solver[6]:d}")
    config.set("Solver", "dx", f"{solver[7]:.6f}")
    config.set("Solver", "boundaryLayers", f"{solver[8]:d}")
    config.set("Solver", "density", f"{solver[9]:.6f}")

    config.add_section("Initialisation")
    config.set("Initialisation", "filling", f"{fill:.6f}")
    config.set("Initialisation", "CoRx", f"{CoR[0]:.6f}")
    config.set("Initialisation", "CoRy", f"{CoR[1]:.6f}")
    config.set("Initialisation", "CoRz", f"{CoR[2]:.6f}")

    config.set("Initialisation", "xi", f"{transAmp[0]:.6f}")
    config.set("Initialisation", "eta", f"{transAmp[1]:.6f}")
    config.set("Initialisation", "zeta", f"{transAmp[2]:.6f}")
    config.set("Initialisation", "theta", f"{rotAmp[0]:.6f}")
    config.set("Initialisation", "psi", f"{rotAmp[1]:.6f}")
    config.set("Initialisation", "chi", f"{rotAmp[2]:.6f}")

    config.set("Initialisation", "xiPhase", f"{transPhase[0]:.6f}")
    config.set("Initialisation", "etaPhase", f"{transPhase[1]:.6f}")
    config.set("Initialisation", "zetaPhase", f"{transPhase[2]:.6f}")
    config.set("Initialisation", "thetaPhase", f"{rotPhase[0]:.6f}")
    config.set("Initialisation", "psiPhase", f"{rotPhase[1]:.6f}")
    config.set("Initialisation", "chiPhase", f"{rotPhase[2]:.6f}")

    config.set("Initialisation", "xiOmega", f"{transOmega[0]:.6f}")
    config.set("Initialisation", "etaOmega", f"{transOmega[1]:.6f}")
    config.set("Initialisation", "zetaOmega", f"{transOmega[2]:.6f}")
    config.set("Initialisation", "thetaOmega", f"{rotOmega[0]:.6f}")
    config.set("Initialisation", "psiOmega", f"{rotOmega[1]:.6f}")
    config.set("Initialisation", "chiOmega", f"{rotOmega[2]:.6f}")

    with open(path, "w") as config_file:
        config.write(config_file)


def readConfig(path):
    import configparser
    config = configparser.ConfigParser()
    config.read(path)

    corx = config.getfloat("Initialisation", "corx")
    cory = config.getfloat("Initialisation", "cory")
    corz = config.getfloat("Initialisation", "corz")

    xi = config.getfloat("Initialisation", "ksi")
    eta = config.getfloat("Initialisation", "eta")
    zeta = config.getfloat("Initialisation", "zeta")
    theta = config.getfloat("Initialisation", "theta")
    psi = config.getfloat("Initialisation", "psi")
    chi = config.getfloat("Initialisation", "xi")

    xiOmega = config.getfloat("Initialisation", "ksiOmega")
    etaOmega = config.getfloat("Initialisation", "etaOmega")
    zetaOmega = config.getfloat("Initialisation", "zetaOmega")
    thetaOmega = config.getfloat("Initialisation", "thetaOmega")
    psiOmega = config.getfloat("Initialisation", "psiOmega")
    chiOmega = config.getfloat("Initialisation", "xiOmega")

    xiPhase = config.getfloat("Initialisation", "ksiPhase")
    etaPhase = config.getfloat("Initialisation", "etaPhase")
    zetaPhase = config.getfloat("Initialisation", "zetaPhase")
    thetaPhase = config.getfloat("Initialisation", "thetaPhase")
    psiPhase = config.getfloat("Initialisation", "psiPhase")
    chiPhase = config.getfloat("Initialisation", "xiPhase")

    #print(config.items("Initialisation"))

    origin = [corx, cory, corz]

    transAmp = [xi, eta, zeta]
    transOmega = [xiOmega, etaOmega, zetaOmega]
    transPhase = [xiPhase, etaPhase, zetaPhase]

    rotAmp = [theta, psi, chi]
    rotOmega = [thetaOmega, psiOmega, chiOmega]
    rotPhase = [thetaPhase, psiPhase, chiPhase]

    return origin, transAmp, transOmega, transPhase, rotAmp, rotOmega, rotPhase