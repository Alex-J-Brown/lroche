import nox


@nox.session(python=["3.9", "3.13"])
def test(session):
    # Install dependencies
    session.install(".[tests]")
    # Run tests
    session.run("pytest")
