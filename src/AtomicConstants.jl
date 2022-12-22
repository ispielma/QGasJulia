"""
Sources
	[1] D. S. Petrov, Ph.D. Thesis, (2003)
	[2] B. Marcelis et al., PRB 70 012701 (2004)
	[3] Le Luo, Ph.D. Thesis, Duke (2009)
	[4] Robert Sylvester Williamson III, Ph.D. Thesis, Wisconson (1997)
"""
module AtomicConstants

    import TOML

    # ======================== GENERAL CONSTANTS ========================

    μ_0 = pi * 4e-7; # magnetic constant
    ϵ_0 = 8.854187817620389e-12; # electric constant
    k_B = 1.3806503e-23; # electric constant
    ħ = 1.0545717253362894e-34; # hbar,
    c = 2.99792458e8; # Speed of light in vacuum
    AMU = 1.660538921e-27; # Atomic Mass Unit

    μ_B = 9.27400968e-24; # Bohr Magneton J/Tesla
    μ_B_Hz_Gauss = 9.27400968e-28 / (2*π*ħ); # Bohr Magneton Hz/Gauss

    a_0 = 5.2917721092e-11; # Bohr Radius
    E_h = 4.35974434e-18; # Hartree energy
    g_earth = 9.80665;  # Gravatational acceleration at the Earths surface

    # ==================== Atomic Transitions Info =====================

    Base.@kwdef mutable struct State # E.g. for a specific term e.g. 5S_{1/2}
        A::Float64=NaN # Hyperfine A coefficient
        B::Float64=NaN # Hyperfine B coefficient

        l::Float64=NaN # orbital angular momentum
        s::Float64=NaN # spin angular momentum
        j::Float64=NaN # Total angular momentum

        g_J::Float64=NaN # Lande g-factor
    end

    Base.@kwdef mutable struct Transition # Properties of a specific transitions
        name::String="" # Nickname of the transition like "D2"
        τ::Float64=NaN # Lifetime
        Γ::Float64=NaN # linewidth
        λ::Float64=NaN # Wavelength
        k::Float64=NaN # wavevector
        ω::Float64=NaN # frequency
    end

    Base.@kwdef struct Atom
        m::Float64=NaN # Atomic Mass
        i::Float64=NaN # Nuclear magnetic moment
        
        # Lande g-factors
        g_S::Float64=NaN
        g_L::Float64=NaN
        g_I::Float64=NaN

        a_s::Float64=NaN # background scattering length
        c_6::Float64=NaN # C6 coefficient * * E_h * a_0**6;

        # Bulk thermodynamics stuff
        T_Melt::Float64=NaN
        T_Boil::Float64=NaN
        
        Vapor_A_Solid::Float64=NaN
        Vapor_A_Liquid::Float64=NaN
        
        Vapor_B_Solid::Float64=NaN
        Vapor_B_Liquid::Float64=NaN

        # states
        states::Dict{String, State} = Dict{String,State}()

        # I don't know how I want to be labeling the transitions
        # But I am sure that this will not be a Time-critical operation
        # So I'll assume "5p_3/2->5s_1/2" will be the format of the keys
        transitions::Dict{String, Transition} = Dict{String,Transition}()
    end

    """
    Single Photon Recoil energy
    """
    Er(at, λ) = ħ^2 * (2*π / λ^2) / (2*at.m)

    """
    Single Photon Recoil energy in Hz
    """
    Er_Hz(at, λ) = Er(at, λ)/ (2*π * ħ)

    """
    Vapor pressure in mBar for a given temperature
    """
    function VaporPressure(at, Temp)
        
        A = Temp > at.T_Melt ? at.Vapor_A_Liquid : at.Vapor_A_Solid
        B = Temp > at.T_Melt ? at.Vapor_B_Liquid : at.Vapor_B_Solid

        10^(A - B/Temp)
    end

    """
    thermal de Broglie wavelength
    """
    Thermal_dB_Wavelength(at, T) = ħ*sqrt(2*npi / (at.m * k_B * T))

    # Now load and process the files in DataFiles/Atoms

    """
    This function converts the TOML data into the atom struct
    """
    function _ProcessAtomTOML(atomdict)
        # Process States
        states = Dict{String,State}()
        statedict = get(atomdict, "states", Dict())

        for (key, value) in statedict
            state = State(
                A=get(value, "A", 0),
                B=get(value, "B", 0),
                l=get(value, "l", 0),
                s=get(value, "s", 0),
                j=get(value, "j", 0),
                g_J=get(value, "g_J", 0)
            )
            states[key] = state
        end

        # Process Transitions
        transitions=Dict{String,Transition}()
        transitiondict = get(atomdict, "transitions", Dict())

        for (key, value) in transitiondict
            transition = Transition(
                name=get(value, "name", ""),
                τ=get(value, "τ", NaN),
                λ=get(value, "λ", NaN),
            )
            transition.Γ = 1 / transition.τ
            transition.k = 2 * π / transition.λ
            transition.ω = c * transition.k

            transitions[key] = transition
        end

        # Create Atom
        atom = Atom(
            m=get(atomdict, "m", NaN)*AMU,
            g_S=get(atomdict, "g_S", NaN),
            g_L=get(atomdict, "g_L", NaN),
            g_I=get(atomdict, "g_I", NaN),
            i=get(atomdict, "i", NaN),
            a_s=get(atomdict, "a_s", NaN)*a_0,
            c_6=get(atomdict, "c_6", NaN),
            T_Melt=get(atomdict, "T_Melt", NaN),
            T_Boil=get(atomdict, "T_Boil", NaN),
            Vapor_A_Solid=get(atomdict, "Vapor_A_Solid", NaN),
            Vapor_A_Liquid=get(atomdict, "Vapor_A_Liquid", NaN),
            Vapor_B_Solid=get(atomdict, "Vapor_A_Solid", NaN),
            Vapor_B_Liquid=get(atomdict, "Vapor_A_Solid", NaN),
            states=states,
            transitions=transitions
        )
    end

    let _sourcepath, _files, _atomnames, _atoms
        _sourcepath = (@__DIR__)  * "/DataFiles/Atoms"

        _files = [file for file in readdir(_sourcepath) if endswith(file, ".toml")]

        _atomnames = Tuple([split(file, ".toml")[1] for file in _files])
        _atoms = [TOML.parsefile(_sourcepath*"/"*file) for file in _files]

        for (name, atom) in zip(_atomnames, _atoms)
            eval(Meta.parse("global $name = _ProcessAtomTOML($atom)"))
        end
    end

end