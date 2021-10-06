//
// Created by Kevin Gori on 02/10/2021.
//

#pragma once

namespace strom {

class XStrom : public std::exception {
public:
    XStrom() noexcept = default;

    explicit XStrom(const std::string s) noexcept : _msg() { _msg = s; }

    XStrom(const XStrom &) noexcept = default;

    virtual ~XStrom() noexcept = default;

    [[nodiscard]] const char *what() const noexcept { return _msg.c_str(); }

private:
    std::string _msg;
};

}// namespace strom