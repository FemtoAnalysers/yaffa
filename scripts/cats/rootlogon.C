void rootlogon() {
    const auto CATS_PATH = gSystem->Getenv("CATS");

    if (!CATS_PATH) {
        std::cout << "\033[33m[WARNING] The CATS environment variable is undefined. It should point to the CATS installation\033[0m" << std::endl;
    } else {
        std::cout << "\033[34m[INFO] Using CATS in " << CATS_PATH << " as per environment variable\033[0m" << std::endl;
        gInterpreter->AddIncludePath(Form("%s/install/include", gSystem->Getenv("CATS")));
    }
}
