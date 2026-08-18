import nox


@nox.session(python=["3.10", "3.11", "3.12"])
def tests(session):
    # Install dependencies
    session.install(".[tests]")
    # Run tests
    session.run("pytest")
