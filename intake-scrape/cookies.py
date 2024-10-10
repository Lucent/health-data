import browser_cookie3

# Load Chrome cookies
cookies = browser_cookie3.firefox(domain_name='.myfitnesspal.com')

# Print each cookie's domain, name, and value
for cookie in cookies:
    print(f"{cookie.domain} | {cookie.name} | {cookie.value}")

